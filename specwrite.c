
#if HAVE_CONFIG_H
#  include <config.h>
#endif

/* System includes */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

/* Starlink includes */
#include "sae_par.h"
#include "star/hds.h"
#include "ndf.h"
#include "ast.h"
#include "ems.h"
#include "mers.h"
#include "prm_par.h"
#include "star/mem.h"

/* Local includes */
#include "specwrite.h"

/* Enable memory cache */
#define USE_MEMORY_CACHE 1

/* Debug prints */

#define SPW_DEBUG 1

/* Largest file name allowed (including path) */
#define MAXFILE 1024

/* Define the number of extensions we support */
#define NEXTENSIONS 34

/* Maximum number of subsystems we can handle 
   We know that ACSIS can have at most 4 spectral subsystems
   and they will not change during a single observation. */

#define MAXSUBSYS 4
static const unsigned int maxsubsys = MAXSUBSYS;

/* Global state variables */

/* This struct gives an overview of the observation state itself */
typedef struct obsData {
  /* number of receptors in this observation */
  unsigned int nrecep;
  size_t receplen;       /* longest receptor name */
  char * recep_name_buff;  /* buffer for receptor names */
  char datadir[MAXFILE+1];   /* current data directory */
  /* observation number, ut date and subscan number */
  unsigned int obsnum;
  unsigned int yyyymmdd;
  unsigned int nsubsys; /* Actual number of subsystems in use */
} obsData;

/* actual data (either in memory or mapped from disk) */
typedef struct specData {
  float  * spectra;  /* Array of data to receive spectra ( nchans x nrecep x nseq) */
  double * receppos; /* Receptor positions (2 x  nrecep x nseq) */
  void   * jcmtstate[NEXTENSIONS];  /* Pointers to JCMTSTATE information */
  float  * bad;      /* array of bad values to easily initialise new time slice */
} specData;

/* This struct contains file information (ndf identifiers, hds locators) */
typedef struct fileInfo {
  unsigned int subscan;  /* current subscan number */
  int indf;   /* NDF identifier for this file */
  HDSLoc * extloc;  /* JCMSTATE extension locator */
  HDSLoc * acsisloc; /* ACSIS extension locator */
  HDSLoc * extlocators[NEXTENSIONS]; /* Locators to each JCMTSTATE component */
  HDSLoc * receppos_loc;   /* Locators to each .MORE.ACSIS.RECEPPOS */  
  int extmapped;  /* JCMTSTATE extension is mapped */
  int acsismapped; /* ACSIS extension is mapped */
} fileInfo;


/* This struct contains all the information required for a single subsystem */
typedef struct subSystem {
  unsigned int index;   /* The subsystem "index" for this subsystem 0, 1, 2, or 3 */
  unsigned int maxsize; /* Maximum allowed size */
  unsigned int cursize; /* Number of sequence steps available */
  unsigned int curseq;  /* current RTS sequence number */
  unsigned int curpos;  /* current index position in the cube (water mark level) */
  unsigned int nchans;  /* number of spectral channels in this subsystem */
  specData     tdata;   /* Actual data */        
  fileInfo     file;    /* File data (can be none if file.infd == NDF__NOID) */
} subSystem;

/* Now pre allocate the required number of subsystems */
subSystem SUBSYS[MAXSUBSYS]; /* subsystem information that is being added to */
subSystem FILEDATA[MAXSUBSYS];
obsData   OBSINFO;           /* Global observation information */

/* Status flags */
unsigned int INPROGRESS = 0; /* true if an observation is in progress */
unsigned int CALLED = 0;     /* has openTS been called at least once */

/* internal prototypes */
static void writeFitsChan( int indf, const AstFitsChan * fitschan, int * status );
static char * getFileName( const char * dir, unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status );
static char * getDirName( const char * dir, unsigned int yyyymmdd, 
			   unsigned int obsnum, int * status );
static char * getFileRoot( unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status );

static char * createSubScanDir( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
				int *status);

static void
openNDF( const obsData * obsinfo, const subSystem * template, subSystem * file,
	 unsigned int nseq, int * status );

static void
closeNDF( subSystem * subsys, int * status );

static void
resizeNDF( const obsData * obsinfo, subSystem * subsys, unsigned int newsize, int * status );

static void
createExtensions( subSystem * subsys, unsigned int size, int * status );

static void
createACSISExtensions( const obsData * obsinfo, subSystem * subsys, unsigned int size, int * status );

static void resizeExtensions( subSystem * subsys, unsigned int newsize, int remap, 
			      int * status );
static void resizeACSISExtensions( subSystem * subsys, unsigned int newsize, int remap, 
			      int * status );
static void closeExtensions( subSystem * subsys, int * status );
static void closeACSISExtensions( subSystem * subsys, int * status );

static void writeRecord( void * basepntr[], unsigned int tindex,
			 const ACSISRtsState * record,
			 int * status );
static void writeRecepPos( const obsData * obsinfo, double * posdata, unsigned int tindex,
			   const ACSISRtsState * record, int * status );
static unsigned int calcOffset( unsigned int nchans, unsigned int maxreceps, unsigned int nrecep, 
				unsigned int tindex, int *status );

static void allocHeaders( subSystem *subsys, unsigned int size, int * status );

static void
resizeResources( const obsData * obsinfo, subSystem * subsys, unsigned int newsize, int * status );
static void
allocResources( const obsData * obsinfo, subSystem *subsys, unsigned int nseq,
		int * status );

static void freeResources ( obsData * obsinfo, subSystem * subsys, int * status);
static void allocPosData( const obsData *obsinfo, subSystem *subsys, unsigned int nseq, int * status );
static void
flushResources( const obsData * obsinfo, subSystem * subsys, int * status );

static void copyCache( const obsData * obsinfo, const subSystem * input,
		       subSystem * output, unsigned int nseq, int * status);

static size_t sizeofHDSType( const char * type, int * status );

static void myRealloc( void **pntr, size_t nbytes, int * status );

#if SPW_DEBUG
static double duration ( struct timeval * tp1, struct timeval * tp2 );
#endif


/* Function to put quotes around a symbol so that we can do
   CPP string concatenation */
#define myxstr(s) mystr(s)
#define mystr(s) #s
#define CHARTYP(s) "_CHAR*" myxstr(s)

/* Threshold for reporting timing anomalise */
#define LONGTIME 2.0

/* Macro to time an event */
#if SPW_DEBUG
#define TIMEME(label,func) { struct timeval tp1; struct timeval tp2; double tpdiff; \
    gettimeofday( &tp1, NULL ); \
    func;			\
    gettimeofday( &tp2, NULL ); \
    tpdiff = duration( &tp1, &tp2 ); \
    if (tpdiff > LONGTIME) printf( ">>>>>>>>>>>>>>>>>" label " took %.3f seconds <<<<<<<<<<<<<<<\n", tpdiff); \
  }
#else
#define TIMEME(func)  func
#endif

/* Macro to check subsys range */
#define CHECKSUBSYS( sub, st ) if ( st == SAI__OK && sub > (MAXSUBSYS -1 )) { \
				      *st = SAI__ERROR; \
				      emsSetu( "SB", sub ); \
				      emsSetu( "MX", (MAXSUBSYS-1) ); \
				      emsRep( " ", "Subsystem out of range. ^SB > ^MX", st ); }

/* Number of dimensions in output NDF */
#define NDIMS 3

/* definitions of dimensions */
#define CHANDIM 0
#define RECDIM  1
#define TDIM    2

/* Define indices for array of mapped pointers to extensions */
#define POL_ANG      0
#define RTS_NUM      1
#define RTS_END      2
#define RTS_TASKS    3
#define SMU_AZ_JIG_X 4
#define SMU_AZ_JIG_Y 5
#define SMU_X        6
#define SMU_Y        7
#define SMU_Z        8 
#define SMU_TR_JIG_X 9
#define SMU_TR_JIG_Y 10
#define TCS_AIRMASS  11
#define TCS_AZ_ANG   12
#define TCS_AZ_AC1   13
#define TCS_AZ_AC2   14
#define TCS_AZ_DC1   15
#define TCS_AZ_DC2   16
#define TCS_AZ_BC1   17
#define TCS_AZ_BC2   18
#define TCS_INDEX    19
#define TCS_SOURCE   20
#define TCS_TR_SYS   21
#define TCS_TR_ANG   22
#define TCS_TR_AC1   23
#define TCS_TR_AC2   24
#define TCS_TR_DC1   25
#define TCS_TR_DC2   26
#define TCS_TR_BC1   27
#define TCS_TR_BC2   28
#define ENVIRO_REL_HUM   29
#define ENVIRO_PRESSURE  30
#define ENVIRO_AIR_TEMP  31
#define ACS_SOURCE_RO    32
#define ACS_DRCONTROL    33

/* Definitions of HDS types associated with ACSISRtsStates struct. All these
   will be created in the file. */
static const char * hdsRecordNames[NEXTENSIONS][2] = 
  {
   { "_DOUBLE", "POL_ANG" },
   { "_INTEGER", "RTS_NUM" },
   { "_DOUBLE", "RTS_END" },
   { CHARTYP(SIZEOF_RTS_TASKS), "RTS_TASKS" },
   { "_DOUBLE", "SMU_AZ_JIG_X" },
   { "_DOUBLE", "SMU_AZ_JIG_Y" },
   { "_DOUBLE", "SMU_X" },
   { "_DOUBLE", "SMU_Y" },
   { "_DOUBLE", "SMU_Z" },
   { "_DOUBLE", "SMU_TR_JIG_X" },
   { "_DOUBLE", "SMU_TR_JIG_Y" },
   { "_DOUBLE", "TCS_AIRMASS" },
   { "_DOUBLE", "TCS_AZ_ANG" },
   { "_DOUBLE", "TCS_AZ_AC1" },
   { "_DOUBLE", "TCS_AZ_AC2" },
   { "_DOUBLE", "TCS_AZ_DC1" },
   { "_DOUBLE", "TCS_AZ_DC2" },
   { "_DOUBLE", "TCS_AZ_BC1" },
   { "_DOUBLE", "TCS_AZ_BC2" },
   { "_INTEGER", "TCS_INDEX" },
   { CHARTYP(SIZEOF_TCS_SOURCE), "TCS_SOURCE" },
   { CHARTYP(SIZEOF_TCS_TR_SYS), "TCS_TR_SYS" },
   { "_DOUBLE", "TCS_TR_ANG" },
   { "_DOUBLE", "TCS_TR_AC1" },
   { "_DOUBLE", "TCS_TR_AC2" },
   { "_DOUBLE", "TCS_TR_DC1" },
   { "_DOUBLE", "TCS_TR_DC2" },
   { "_DOUBLE", "TCS_TR_BC1" },
   { "_DOUBLE", "TCS_TR_BC2" },
   { "_REAL", "ENVIRO_REL_HUM" },
   { "_REAL", "ENVIRO_PRESSURE" },
   { "_REAL", "ENVIRO_AIR_TEMP" },
   { CHARTYP(SIZEOF_ACS_SOURCE_RO), "ACS_SOURCE_RO" },
   { "_INTEGER", "ACS_DRCONTROL" },
  };

/* Somewhere to store the precomputed sizes of each HDS element */
static size_t hdsRecordSizes[NEXTENSIONS];

/* Extension support */

/* Name of STATE and ACSIS extensions - some could almost be shared with SCUBA2... */
#define STATEEXT   "JCMTSTATE"
#define STATEEXTTYPE "RTS_ARR"
#define ACSISEXT   "ACSIS"
#define ACSISEXTTYP "ACSIS_COMP"

/*********************** NDF "cube" FILE *************************************/

/* Number of sequences to increment file size by if it runs out of space */

/* May want this parameter to be settable as the number of spectra that
   we expect in a given time period */
#define MAXRECEP   16
#define MAXRATE    20
#define PRESIZETIME 10
#define NGROW  (MAXRATE * PRESIZETIME)

/* Number of bytes we should write before opening a new file */
#define MAXBYTES ( 500 * 1024 * 1024 )

/* maximum number of sequence steps we can get out of sequence */
#define MAXSEQERR  5


/*
*+
*  Name:
*     acsSpecOpenTS

*  Purpose:
*     Open NDF file for writing time series spectra

*  Invocation:
*     acsSpecOpenTS( const char * dir, unsigned int yyyymmdd, 
*                    unsigned int obsnum, unsigned int nrecep,
*                    unsigned int nsubsys, const size_t nchans[], 
*                    unsigned int nseq, const char *recepnames[],
*                    int * status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to prepare the file for output. It must
*     be called before calling acsSpecWriteTS. It will be pre-sized to receive
*     spectra.

*  Arguments:
*     dir = const char * (Given)
*        Directory to write the file.
*     yyyymmdd = unsigned int (Given)
*        UT date in YYYYMMDD format. Used to construct the file name.
*     obsnum = unsigned int (Given)
*        Current observation number. Used to construct the file name.
*     nrecep = unsigned int (Given)
*        Number of receptors participating in this observation.
*     nsubsys = unsigned int (Given)
*        Number of subsystems to be written in this observation.
*     nchans[] = size_t (Given)
*        Number of channels in each subsystem. Should have at least "nsubsys"
*        elements.
*     nseq = unsigned int (Given)
*        Initial number of sequence steps to use to presize the data file.
*        If 0 is given, an internal default will be used.
*     recepnames[] = const char*[] (Given)
*        Names of each receptor in feed order.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.
*     05-APR-2006 (TIMJ):
*        Use structured globals

*  Notes:
*     - Currently only one spectrum file can be open for write at any
*     given time. It is an error for this function to be called whilst
*     a file is open. Call acsSpecClose to close the file.
*     - The file created by this routine will be of the form
*       aYYYYMMDD_NNNNN.sdf where NNNNN is the zero padded observation number.
*     - The resulting file will be an HDS container containing an NDF
*     for each subsystem. These NDFs will be dimensioned as nchans*time.

*  Copyright:
*     Copyright (C) 2006 Particle Physics and Astronomy Research Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*
*     This program is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public
*     License along with this program; if not, write to the Free
*     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*     MA 02111-1307, USA

*  Global Variables:
*     OBSINFO
*     SUBSYS
*     hdsRecordName      (readonly)
*     NEXTENSIONS        (readonly)

*-
*/

void
acsSpecOpenTS( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
	       unsigned int nrecep, unsigned int nsubsys, 
	       const size_t nchans[], unsigned int nseq,
	       const char *recepnames[], int * status ) {

  char * cpos = NULL;          /* offset into string */
  unsigned int i;              /* Loop counter */
  unsigned int j;              /* Loop counter */
  char * sdir = NULL;          /* subscan directory */
  size_t len;                  /* temp length */
  size_t receplen;             /* Length of longest receptor name */
  unsigned int nperseq;        /* Number of elements per sequence */
  char * recep_name_buff = NULL; /* Malloced buffer for receptor names */
  subSystem * subsys = NULL;   /* Individual subsystem */

  /* Return immediately if status is bad */
  if (*status != SAI__OK) return;

  /* Validate the subsystem count */
  if (nsubsys > maxsubsys) {
    *status = SAI__ERROR;
    emsSetu("NIN", nsubsys);
    emsSetu("MAX", maxsubsys);
    emsRep("HDS_SPEC_OPENTS_ERR0",
	   "acsSpecOpenTS: number of subsystems supplied (^NIN) exceeds expected maximum of ^MAX", status);
    return;
  }

  /* Check to see if we've already been called */
  if (INPROGRESS != 0) {
    *status = SAI__ERROR;
    emsRep("HDS_SPEC_OPENTS_ERR1",
	   "acsSpecOpenTS called, yet an observation is already in progress", status);
    return;
  }

  /* Pre-calculate the byte sizes of the extensions */
  /* Loop and create  */
  if (!CALLED) {
    CALLED = 1;
    for (i=0; i < NEXTENSIONS; i++ ) {
      /* work out the number of bytes per element */
      hdsRecordSizes[i] = sizeofHDSType( hdsRecordNames[i][0], status );
      if (*status != SAI__OK) break;
    }
  }

  /* Populate the observation info structure */
  OBSINFO.yyyymmdd = yyyymmdd;
  OBSINFO.obsnum   = obsnum;
  OBSINFO.nsubsys  = nsubsys;
  OBSINFO.nrecep  = nrecep;

  /* Create the directory to receive each subscan */
  sdir = createSubScanDir( dir, yyyymmdd, obsnum, status );

  /* Store the resulting directory */
  strncpy(OBSINFO.datadir, sdir, MAXFILE);
  (OBSINFO.datadir)[MAXFILE] = '\0';


  /* we need to store the receptor names somewhere locally */
  /* just take the simple approach of a buffer of strings of equal
     length that we can datPutC */
  /* find the longest receptor name */
  receplen = 0;
  if (recepnames != NULL) {
    for (i=0; i < nrecep; i++) {
      len = strlen( recepnames[i] );
      if (receplen < len ) receplen = len;
    }
  }

  /* Allocate memory for the receptor names */
  if (receplen > 0 ) {
    recep_name_buff = starMalloc( receplen * nrecep );

    /* copy characters into buffer */
    if (recep_name_buff != NULL) {
      cpos = recep_name_buff;
      for (i=0; i<nrecep; i++) {
	cnfExprt( recepnames[i], cpos, receplen);
	cpos += receplen;
      }
    }
  }

  /* Store receptor information */
  OBSINFO.receplen = receplen;
  OBSINFO.recep_name_buff = recep_name_buff;

  /* Need an NDF per subsystem */
  for (i = 0; i < nsubsys; i++) {

    /* Select the subsystem to modify */
    subsys = &(SUBSYS[i]);

    /* zero it out */
    memset( subsys, 0, sizeof(*subsys));

    /* Some intialisation */
    subsys->nchans = nchans[i];
    subsys->index = i;

    /* Number of data values per sequence */
    nperseq = nrecep * nchans[i];

    /* Calculate the number of sequence steps that we are allowed to grow
       before opening new file. */
    subsys->maxsize = MAXBYTES / ( nperseq * SIZEOF_FLOAT );

    /* Allocate a cache of bad values to simplify initialisation */
    subsys->tdata.bad = starMalloc( nperseq * SIZEOF_FLOAT );
    if (subsys->tdata.bad == NULL) {
      if (*status == SAI__OK) {
	*status = SAI__ERROR;
	emsRep(" ","Unable to allocate memory for bad value cache", status );
	break;
      }
    }
    for (j=0; j<nperseq; j++) {
      (subsys->tdata.bad)[j] = VAL__BADR;
    }

    /* Allocate resources for this subsystem */
    allocResources( &OBSINFO, subsys, subsys->maxsize, status );

  }

  /* indicate that an observation is now ready */
  INPROGRESS = 1;

  /* Complete */
  return;

}

/*
*+
*  Name:
*     acsSpecWriteTS

*  Purpose:
*     Write a spectrum to an HDS time-series file.

*  Invocation:
*     acsSpecWriteTS( unsigned int subsys, const float spectrum[], 
*                     const ACSISRtsState * record, const AstFitsChan * freq,
*                     int *status);

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to write a spectrum to an HDS container
*     that had previously been opened using acsSpecOpen.

*  Arguments:
*     subsys = unsigned int (Given)
*        Subsystem used for this spectrum (start counting at 0).
*        Much match the nchans[] array given to acsSpecOpenTS.
*     spectrum = float[nchan] (Given)
*        Spectrum itself.
*     record = const ACSISRtsState * (Given)
*        Header information associated with this spectrum.
*     freq = const AstFitsChan * (Given)
*        Spectral coordinate information for this subsystem.
*        If non-NULL, the spectral information is extracted
*        and stored internally for this subsystem. Should be
*        supplied with the first spectrum from each subsystem.
*        Subsequent calls can pass in NULL and the cached
*        world coordinates will be used.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.
*     05-APR-2006 (TIMJ):
*        Use structured globals

*  Notes:
*     - Must have previously called acsSpecOpenTS.

*  Copyright:
*     Copyright (C) 2006 Particle Physics and Astronomy Research Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*
*     This program is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public
*     License along with this program; if not, write to the Free
*     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*     MA 02111-1307, USA

*-
*/

void
acsSpecWriteTS( unsigned int subsysnum, const float spectrum[], 
		const ACSISRtsState* record,
		const AstFitsChan * freq, int * status ) {

  float * data; /* local copy of mapped pointer to spectrum */
  unsigned int offset;         /* offset into data array */
  int seqinc = 0;              /* did we increment sequence number? */
  unsigned int ngrow;          /* number of time slices to grow file */
  unsigned int reqnum;         /* number of time slices indicates by RTS sequence */
  AstSpecFrame * template;     /* Template of a specFrame */
  AstFrameSet  * foundframe;   /* Frameset from template to found frame */
  unsigned int tindex;         /* Position in sequence array for this sequence */
  unsigned int *rtsseqs;       /* Pointer to array of sequence numbers */
  unsigned int max_behind;     /* How many sequence steps to look behind */
  int found = 0;               /* Did we find the sequence? */
  unsigned int i;              /* loop counter */
  double * posdata;
  void ** recdata;
  unsigned int startind;
  subSystem * subsys;

  if (*status != SAI__OK) return;

  /* make sure that the subsys number is in range */
  if ( subsysnum >= maxsubsys ) {
    *status = SAI__ERROR;
    emsSetu("IN", subsysnum);
    emsSeti("MAX", maxsubsys-1);
    emsRep(" ","acsSpecWriteTS: Supplied subsystem number (^IN) exceeds max allowed (^MAX)", status);
    return;
  }

  /* Check to see if we've already been called */
  if (!INPROGRESS) {
    *status = SAI__ERROR;
    emsRep("HDS_SPEC_WRITETS_ERR1",
	   "acsSpecWriteTS called, yet an observation has not been initialised", status);
    return;
  }

  /* Check feed range */
  if ( record->acs_feed >= OBSINFO.nrecep ) {
    *status = SAI__ERROR;
    emsSetu( "NR", OBSINFO.nrecep - 1 );
    emsSetu( "FEED", record->acs_feed );
    emsRep( " ", "acsSpecWriteTS called, yet the feed number (^FEED) exceeds the expected number (^NR)",
	    status );
  }

  /* Get local copy of subsystem from global */
  subsys = &(SUBSYS[subsysnum]);

  /* first thing to do is determine whether this sequence number is new or old */

  if (subsys->curpos == 0) {
    /* have not written anything yet so the correct place to write this
       sequence is at index 0 in the data array (assuming cursize[] is
       large enough) */
    seqinc = 1;
    tindex = 0;
    printf("First spectrum for this file for this subsystem\n");

  } else {
    /* if the supplied value is the most recent value then we do not 
       need to search */

    if (record->rts_num == subsys->curseq) {

      tindex = subsys->curpos - 1;
      /* printf("Reusing sequence %u at index %u\n", curseq[subsys], tindex);*/

    } else {

      /* pointer to array of rts sequence numbers */
      rtsseqs = (subsys->tdata.jcmtstate)[RTS_NUM];

      /* search back through the list */
      /* For efficiency, assume that we can't be more than a certain
	 number of sequence steps behind. */
      max_behind = ( subsys->curpos < MAXSEQERR ? subsys->curpos : MAXSEQERR );
      startind = subsys->curpos;
      found = 0;
      for (i = 1; i < max_behind; i++) {
	if (record->rts_num == rtsseqs[startind-i]) {
	  /* found the sequence */
	  tindex = startind - i;
	  found = 1;
	  break;
	}
      }

      if (found == 0) {
	printf("Did not find sequence %u\n", record->rts_num);
	/* did not find this sequence number so it is a new one */
	tindex = subsys->curpos; /* curpos will be incremented */
	seqinc = 1;
      } else {
	/* going back in time so this may not be efficient */
	printf("Sequence %u matches index %u. Previous seq num=%u\n", 
	       record->rts_num, tindex, subsys->curseq);
	subsys->curseq = record->rts_num;

      }
    }
  }
  

  /* if the sequence number has incremented we need to increase the t-axis counter */

  if (seqinc) {

    /* store the new value */
    subsys->curseq = record->rts_num;

    /* increment the counters value */
    (subsys->curpos)++;

    /* See if we need to grow */
    /* We can grow either because we have suddenly realised we don't fit *or* because
       we have been told how many sequence steps to expect - calculate the required
       number to extend. (but we know at least 1) */
    reqnum = 1;

    if ( record->rts_endnum > record->rts_num ) {
      reqnum = record->rts_endnum - record->rts_num + 1;
    }

    if ( (subsys->curpos + reqnum - 1) > subsys->cursize ) {
      /* resize NDF and all data arrays */

      /* work out how much to grow:
	 - use requested size if given and more than 1
	 - else use NGROW
	 - make sure we grow by at least subsys->curpos-cursize[subsys]
      */

      if (reqnum > 1) {
	ngrow = reqnum + subsys->curpos - subsys->cursize - 1;
      } else {
	ngrow = NGROW;
      }

      /* if this will put us over the limit, close the NDF and reopen a new
	 one */
      if ( (subsys->cursize + ngrow ) > subsys->maxsize ) {

	/* this will clear cursize but we need to make sure that it reflects the
	   actual number of elements written so far, not the position we were
	   going to write to.*/
	subsys->curpos--;

	/* flush what we have to disk */
	flushResources( &OBSINFO, subsys, status );

#if USE_MEMORY_CACHE
	/* always want to make sure we allocate the max amount of memory */
	ngrow = subsys->maxsize;
#endif
	allocResources( &OBSINFO, subsys, ngrow, status );

	/* indicate that we are starting at the beginning with the next spectrum */
	tindex = 0;

      } else {
	printf("Cursize: %u Curpos: %u ngrow: %u maxsize: %u\n",
	       subsys->cursize, subsys->curpos, ngrow, subsys->maxsize); 
	/* Resize the NDF */
	resizeResources( &OBSINFO, subsys, ngrow, status );
      }
    }
  }

  /* copy in the data */
  if (*status == SAI__OK) {
    
    /* Calculate offset into array - number of spectra into the array times number of
       channels per spectrum. */
    offset = calcOffset( subsys->nchans, OBSINFO.nrecep,
			 record->acs_feed, tindex, status );

    data = (subsys->tdata.spectra);
    recdata = (subsys->tdata.jcmtstate);
    posdata = (subsys->tdata.receppos);

    memcpy( &(data[offset]), spectrum, subsys->nchans*SIZEOF_FLOAT );

    /* Store record data and receptor positions. Base record only updates each
       sequence step but recepot position should be written for all records. */
    if (seqinc) writeRecord( recdata, tindex, record, status );
    writeRecepPos( &OBSINFO, posdata, tindex, record, status );

  }

  return;
}

/*
*+
*  Name:
*     acsSpecCloseTS

*  Purpose:
*     Write FITS header and close HDS file.

*  Invocation:
*     acsSpecCloseTS( const AstFitsChan * fits[], int *status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to close the file after all spectra have been
*     written. The FITS header is written.

*  Arguments:
*     fits[] = const AstFitsChan * (Given)
*        Array of FITS headers. One per subsystem.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.
*     20-APR-2006 (TIMJ):
*        Use structured globals.

*  Notes:
*     - Must have previously called acsSpecOpenTS.
*     - File is resized to the actual number of spectra written.

*  Copyright:
*     Copyright (C) 2006 Particle Physics and Astronomy Research Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*
*     This program is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public
*     License along with this program; if not, write to the Free
*     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*     MA 02111-1307, USA

*-
*/

void
acsSpecCloseTS( const AstFitsChan * fits[], int * status ) {

  unsigned int i;           /* Loop counter */
  int found = 0;            /* Found an open NDF? */
  subSystem * subsys;       /* specific subsystem */

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Check to see if we've already been called */
  if (!INPROGRESS) {
    *status = SAI__ERROR;
    emsRep("HDS_SPEC_CLOSETS_ERR1",
	   "acsSpecCloseTS called, yet an observation has not been initialised", status);
    return;
  }


  /* Loop over each NDF to write fits header and to close it */
  found = 0;
  for (i = 0; i < maxsubsys; i++) {
    /* Get local copy of subsystem from global */
    subsys = &(SUBSYS[i]);
    if ( subsys->file.indf != NDF__NOID || subsys->tdata.spectra != NULL) {
      found = 1;
      flushResources( &OBSINFO, subsys, status);
    }
    freeResources( &OBSINFO, subsys, status );
  }

  /* report error if not found any open NDFs */
  if (*status == SAI__OK && !found) {
    *status = SAI__ERROR;
    emsRep(" ", "acsSpecCloseTS: Failed to find open NDF components", status );
  }

  /* Now need to open all the files that we have opened previously and adjust
     all the FITS headers. This is going to be problematic for the ones that
     are related to sequence values.
  */

  /* Now need to write out the .ok file with all the files that we have opened */


  /* Force globals to be reset */
  for (i = 0; i < maxsubsys; i++) {
    memset( &(SUBSYS[i]), 0, sizeof(SUBSYS[i]) );
  }

  /* reset progress */
  INPROGRESS = 0;

  if (*status != SAI__OK) {
    emsRep( " ", "Error closing Spectrum file", status );
  }

}

/*********************** HELPER FUNCTIONS (PRIVATE) *****************************/

/* Internal function to write a FitsChan to an NDF - taken from
   kpgPtfts */

static void writeFitsChan( int indf, const AstFitsChan * fits, int * status ) {

  char card[81];            /* A single FITS header card */
  HDSLoc * fitsloc = NULL;  /* Locator to FITS extension */
  char * fpntr;             /* Pointer to mapped FITS header */
  unsigned int i;           /* Loop counter */
  unsigned int ncards;      /* Number of header cards */
  size_t nchars;            /* Actual size of FITS extension */
  int result;               /* Result from astFindFits */
  int extdims[1];

  if (*status != SAI__OK) return;

  /* Find out how many header cards we have */
  ncards = astGetI( fits, "Ncard" );

  /* Rewind the FitsChan */
  astClear( fits, "Card" );
    
  /* Create FITS extension */
  extdims[0] = ncards;
  ndfXnew(indf, "FITS", "_CHAR*80", 1, extdims, &fitsloc, status );

  /* Loop over all cards, inserting into extension */
  datMapV( fitsloc, "_CHAR*80", "WRITE", (void**)&fpntr, &nchars, status );

  if (*status == SAI__OK) {
    if ( ncards != nchars ) {
      *status = SAI__ERROR;
      emsSetu( "DM", nchars );
      emsSetu( "SZ",  ncards );
      emsRep("HDS_SPEC_CLOSE_ERR2",
	     "Bizarre error whereby number of cards in mapped FITS header (^DM) differs from number requested (^SZ)", status );
    }
  }

  if (*status == SAI__OK) {
    for (i = 1; i <= ncards; i++) {
      result = astFindFits( fits, "%f", card, 1 );
      if (result) {
	strncpy( fpntr, card, 80 );
	fpntr += 80;
      } else {
	break;
      }
    }
  }

  /* Cleanup */
  datUnmap( fitsloc, status );
  datAnnul( &fitsloc, status );

  if (*status != SAI__OK)
    emsRep(" ", "Error writing FITS information.", status );
}

/* Form the file name
   - returns a pointer to static memory.
 */

static char * getFileName( const char * dir, unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status ) {

  static char filename[MAXFILE]; /* buffer for filename - will be returned */
  int flen;                        /* Length of string */
  char * root;                   /* Root filename */

  if (*status != SAI__OK) return NULL;
  root = getFileRoot( yyyymmdd, subsys, obsnum, subscan, status );
  if (*status != SAI__OK) return NULL;

  /* Form the file name - assume posix filesystem */
  flen = snprintf(filename, MAXFILE, "%s/%s.sdf", dir, root );

  if (flen >= MAXFILE) {
    *status = SAI__ERROR;
    emsSeti("SZ", MAXFILE );
    emsSetu("N", obsnum );
    emsSetu("UT", yyyymmdd );
    emsSetu("SS", subscan );
    emsRep("HDS_SPEC_OPEN_ERR1",
	   "Error forming filename. Exceeded buffer size of ^SZ chars for scan ^N subscan ^SS on UT ^UT", status );
    return NULL;
  }

  return filename;
}

/* Get the root of the name. No directory and no suffix.
   If subscan is zero, the subscan is not included (useful for building directory name).
   Returns static memory.
 */

static char * getFileRoot( unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status ) {

  static char rootname[MAXFILE]; /* buffer for filename - will be returned */
  int flen = 0;                  /* Length of string from snprintf */

  if (*status != SAI__OK) return NULL;

  /* Check subsys */
  CHECKSUBSYS(subsys, status );

  /* Form the file name - assume posix filesystem */
  if (*status == SAI__OK) {
    if (subscan > 0) {
      flen = snprintf(rootname, MAXFILE, "a%u_%05u_%02u_%04u", yyyymmdd, obsnum, subsys, subscan );
    } else {
      flen = snprintf(rootname, MAXFILE, "a%u_%05u", yyyymmdd, obsnum );
    }
  }

  if (flen >= MAXFILE && *status == SAI__OK) {
    *status = SAI__ERROR;
    emsSeti("SZ", MAXFILE );
    emsSetu("N", obsnum );
    emsSetu("UT", yyyymmdd );
    emsRep("HDS_SPEC_OPEN_ERR1",
	   "Error forming filename. Exceeded buffer size of ^SZ chars for scan ^N on UT ^UT", status );
    return NULL;
  }

  return rootname;
}

/* Form the directory name
   To comply with the SCUBA-2 ICD we only use the observation number when
   constructing the directory name.

   - returns a pointer to static memory.
 */

static char * getDirName( const char * dir, unsigned int yyyymmdd,
			  unsigned int obsnum, int * status ) {

  static char dirname[MAXFILE]; /* buffer for dirname - will be returned */
  int flen;                        /* Length of string */
  char * root;                   /* Root filename */

  if (*status != SAI__OK) return NULL;

  /* Form the file name - assume posix filesystem */
  flen = snprintf(dirname, MAXFILE, "%s/%05u", dir, obsnum );

  if (flen >= MAXFILE) {
    *status = SAI__ERROR;
    emsSeti("SZ", MAXFILE );
    emsSetu("N", obsnum );
    emsSetu("UT", yyyymmdd );
    emsRep("HDS_SPEC_OPEN_ERR1",
	   "Error forming directory name. Exceeded buffer size of ^SZ chars for scan ^N on UT ^UT", status );
    return NULL;
  }

  return dirname;
}


/* Open a directory to contain the subscans. Returns pointer to static
   memory containing the directory name. This will be overwritten in a
   subsequent call. */

static char *
createSubScanDir( const char * rootdir, unsigned int yyyymmdd, unsigned int obsnum,
		  int * status ) {

  char * dir;
  int st;      /* mkdir status */
  mode_t mode;

  if (*status != SAI__OK) return NULL;

  dir = getDirName( rootdir, yyyymmdd, obsnum, status );

  if ( *status == SAI__OK) {
    mode = S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
    st = mkdir( dir, mode );

    if (st == -1) {
      /* already existing is okay (for now) */
      if (errno != EEXIST) {
	*status = SAI__ERROR;
	emsSyser( "MESSAGE", errno );
	emsSetc( "DIR", rootdir );
	emsRep( " ",
		"Error opening subscan directory '^DIR' : ^MESSAGE", status );
      }
    }
  }

  if (*status == SAI__OK) {
    return dir;
  } else {
    return NULL;
  }

}

/* Open an NDF file for that subsystem
   - note that we pass in a subsystem struct for the template
     and populate a subsystem struct corresponding to the file on disk
     The "file" struct should have been initialised since subscan number
     is extracted from it. "template" can be the same struct as "file"
     so you can still use a file as memory cache.
   - The template is simply used to work out how big the output file should
     be.
*/

static void
openNDF( const obsData * obsinfo, const subSystem * template, subSystem * file,
	 unsigned int nseq, int * status ) {

  void *datapntrs[] = { NULL };/* Array of mapped pointers for ndfMap */
  int itemp;                   /* Temp integer */
  int lbnd[NDIMS];             /* Lower pixel bounds */
  char *ndfname;               /* NDF filename */
  unsigned int ngrow;          /* Initial size to grow array */
  int place;                   /* NDF placeholder */
  int ubnd[NDIMS];             /* upper pixel bounds */
  char * history[1] = { "ACSIS Data Acquistion" };


  if (*status != SAI__OK) return;

  if (nseq > template->maxsize) {
    *status = SAI__ERROR;
    emsSetu( "NS", nseq);
    emsSetu( "MAX", template->maxsize );
    emsSetu( "MB", MAXBYTES / (1024*1024));
    emsRep(" ", "openNDF: Unable to open NDF since the number of sequences to be stored (^NS) already exceeds the maximum allowed (^MAX seq equivalent to ^MB megabytes)", status);
    return;
  }

  /* increment subscan number for this file - we assume it has got the correct previous
     value - should this explicitly come from the template? */
  (file->file.subscan)++;

#if SPW_DEBUG
  printf("Opening file associated with subsystem %d (subscan %u)\n",template->index, file->file.subscan);
#endif
  /* Name the NDF component */
  ndfname = getFileName(obsinfo->datadir, obsinfo->yyyymmdd, template->index, obsinfo->obsnum,
			file->file.subscan, status);

  /* Calculate bounds */
  ngrow = (nseq > 0 ? nseq : NGROW );
  lbnd[RECDIM] = 1;
  lbnd[CHANDIM] = 1;
  lbnd[TDIM] = 1;
  ubnd[RECDIM] = obsinfo->nrecep;
  ubnd[CHANDIM] = template->nchans;
  ubnd[TDIM] = ngrow;

  /* Make sure that "file" gets values from "template" */
  file->index = template->index;
  file->maxsize = template->maxsize;
  file->nchans = template->nchans;

#if SPW_DEBUG
  printf("Opening NDF file '%s' to default size of %u sequence steps\n", ndfname, ngrow);
#endif

  /* create the NDF */
  ndfPlace( NULL, ndfname, &place, status );
  ndfNew( "_REAL", NDIMS, lbnd, ubnd, &place, &(file->file.indf), status );

  /* Update the cursize[] array and the nchans array */
  file->cursize = ubnd[TDIM] - lbnd[TDIM] + 1;

  /* History component */
  ndfHcre( file->file.indf, status );
  ndfHput("NORMAL","ACSIS-DA (V" PACKAGE_VERSION ")", 1, 1, history,
	  0, 0, 0, file->file.indf, status );
  ndfHsmod( "DISABLED", file->file.indf, status );

  /* Map the data array */
  ndfMap(file->file.indf, "DATA", "_REAL", "WRITE", datapntrs, &itemp, status );

  /* Store the pointer */
  file->tdata.spectra = datapntrs[0];
  file->curpos = 0;
  file->curseq = 0;

  /* Create the ACSIS extension that contains the receptor names and
     positions */
  createACSISExtensions( obsinfo, file, ngrow, status );

  /* Also need to create the header arrays and map those ! */
  createExtensions( file, ngrow, status );

  return;
}

/* 
   Create the extensions of specified size.
   Pointers stored in extdata.
   HDS Locators stored in extlocators.
*/

static void
createExtensions( subSystem * subsys, unsigned int size, int * status ) {

  int j;
  hdsdim dim[1];
  size_t ndim = 1;

  if (*status != SAI__OK) return;

  if (subsys->file.extmapped) {
    *status = SAI__ERROR;
    emsRep(" ", "createExtensions: Attempting to create extensions that are already mapped", status );
    return;
  }

  /* Initial size */
  dim[0] = size;

  /* Create the extension */
  ndfXnew( subsys->file.indf, STATEEXT, STATEEXTTYPE, 0, NULL, &(subsys->file.extloc), status ); 

  /* Loop and create. Can initialise HDS locator array safely */
  for (j=0; j < NEXTENSIONS; j++ ) {
    (subsys->file.extlocators)[j] = NULL;
    datNew( subsys->file.extloc, hdsRecordNames[j][1], hdsRecordNames[j][0],
	    ndim, dim, status );

    datFind( subsys->file.extloc, hdsRecordNames[j][1], &((subsys->file.extlocators)[j]), status );

    datMap( (subsys->file.extlocators)[j], hdsRecordNames[j][0], "WRITE",
	    ndim, dim, &((subsys->tdata.jcmtstate)[j]), status );
    if ( *status != SAI__OK ) break;

  }

  if (*status == SAI__OK) subsys->file.extmapped = 1;

  if (*status != SAI__OK)
    emsRep(" ", "Error creating JCMT state extension", status );

}

/*
  Resize the extensions to the supplied value.
  If remap is false, the arrays will not be remapped (so call at end to resize
  before annulling locators 
*/

static void
resizeExtensions( subSystem * subsys, unsigned int newsize, 
		  int remap, int * status ) {

  int j;
  hdsdim dim[1];
  size_t ndim = 1;

  if (*status != SAI__OK) return;

  dim[0] = newsize;

  /* Do all the unmapping. Then all the resizing then all the mapping */

  for (j=0; j < NEXTENSIONS; j++ ) {

    datUnmap( subsys->file.extlocators[j], status );
    if ( *status != SAI__OK ) break;
  }

  for (j=0; j < NEXTENSIONS; j++ ) {

    /* resize */
    datAlter( (subsys->file.extlocators)[j], 1, dim, status);
    if ( *status != SAI__OK ) break;
  }

  if (remap) {
    for (j=0; j < NEXTENSIONS; j++ ) {
      
      /* remap - assume this should be done after resizing all */
      datMap( (subsys->file.extlocators)[j], hdsRecordNames[j][0], "WRITE",
	      ndim, dim, &((subsys->tdata.jcmtstate)[j]), status );
      if ( *status != SAI__OK ) break;

    }
  }

  if (*status != SAI__OK)
    emsRep(" ", "Error resizing JCMT state extension", status );

}

/* Close down the extensions and free resources */

static void closeExtensions( subSystem * subsys, int * status ) {

  int j;

  if ( *status != SAI__OK ) return;

  if (subsys->curpos > 0) {
    resizeExtensions( subsys, subsys->curpos, 0, status );
  }

  /* Free locators */
  for (j=0; j < NEXTENSIONS; j++) {
    datAnnul( &((subsys->file.extlocators)[j]), status );
  }

  /* Close extension */
  datAnnul( &(subsys->file.extloc), status );

  /* indicate that we are closed down */
  subsys->file.extmapped = 0;

  /* delete the extension if we never wrote to it */
  if (subsys->curpos == 0) {
    ndfXdel(subsys->file.indf, STATEEXT,status);
  }

  if (*status != SAI__OK)
    emsRep(" ", "Error closing JCMT state extension", status );
}

/* Write ACSISRtsState to file */

static void writeRecord( void * basepntr[], unsigned int frame,
			 const ACSISRtsState * record,
			 int * status ) {
  /* Can not think of anything clever to do */
  if ( *status != SAI__OK ) return;

  /* now copy */
  ((double *)basepntr[POL_ANG])[frame] = record->pol_ang;
  ((int *)basepntr[RTS_NUM])[frame] = record->rts_num;
  ((double *)basepntr[RTS_END])[frame] = record->rts_end;

  cnfExprt( record->rts_tasks,
	   (char *)basepntr[RTS_TASKS]+
	    SIZEOF_RTS_TASKS*frame,
	    SIZEOF_RTS_TASKS );

  ((double *)basepntr[SMU_AZ_JIG_X])[frame] = record->smu_az_jig_x;
  ((double *)basepntr[SMU_AZ_JIG_Y])[frame] = record->smu_az_jig_y;
  ((double *)basepntr[SMU_X])[frame] = record->smu_x;
  ((double *)basepntr[SMU_Y])[frame] = record->smu_y;
  ((double *)basepntr[SMU_Z])[frame] = record->smu_z;
  ((double *)basepntr[SMU_TR_JIG_X])[frame] = record->smu_tr_jig_x;
  ((double *)basepntr[SMU_TR_JIG_Y])[frame] = record->smu_tr_jig_y;
  ((double *)basepntr[TCS_AIRMASS])[frame] = record->tcs_airmass;
  ((double *)basepntr[TCS_AZ_ANG])[frame] = record->tcs_az_ang;
  ((double *)basepntr[TCS_AZ_AC1])[frame] = record->tcs_az_ac1;
  ((double *)basepntr[TCS_AZ_AC2])[frame] = record->tcs_az_ac2;
  ((double *)basepntr[TCS_AZ_DC1])[frame] = record->tcs_az_dc1;
  ((double *)basepntr[TCS_AZ_DC2])[frame] = record->tcs_az_dc2;
  ((double *)basepntr[TCS_AZ_BC1])[frame] = record->tcs_az_bc1;
  ((double *)basepntr[TCS_AZ_BC2])[frame] = record->tcs_az_bc2;
  ((int *)basepntr[TCS_INDEX])[frame] = record->tcs_index;
  cnfExprt ( record->tcs_source,
	     (char *)basepntr[TCS_SOURCE]+SIZEOF_TCS_SOURCE*frame, 
	     SIZEOF_TCS_SOURCE );

  cnfExprt ( record->tcs_tr_sys,
	     (char *)basepntr[TCS_TR_SYS]+SIZEOF_TCS_TR_SYS*frame, 
	     SIZEOF_TCS_TR_SYS );

  ((double *)basepntr[TCS_TR_ANG])[frame] = record->tcs_tr_ang;
  ((double *)basepntr[TCS_TR_AC1])[frame] = record->tcs_tr_ac1;
  ((double *)basepntr[TCS_TR_AC2])[frame] = record->tcs_tr_ac2;
  ((double *)basepntr[TCS_TR_DC1])[frame] = record->tcs_tr_dc1;
  ((double *)basepntr[TCS_TR_DC2])[frame] = record->tcs_tr_dc2;
  ((double *)basepntr[TCS_TR_BC1])[frame] = record->tcs_tr_bc1;
  ((double *)basepntr[TCS_TR_BC2])[frame] = record->tcs_tr_bc2;

  cnfExprt( record->acs_source_ro,
	    (char*)basepntr[ACS_SOURCE_RO]+ SIZEOF_ACS_SOURCE_RO*frame,
	    SIZEOF_ACS_SOURCE_RO );

  ((int *)basepntr[ACS_DRCONTROL])[frame] = record->acs_drcontrol;

  ((float *)basepntr[ENVIRO_AIR_TEMP])[frame] = record->enviro_air_temp;
  ((float *)basepntr[ENVIRO_PRESSURE])[frame] = record->enviro_pressure;
  ((float *)basepntr[ENVIRO_REL_HUM])[frame] = record->enviro_rel_hum;

}

/* Create the .MORE.ACSIS extensions (that are not JCMT state structure members) */


static void
createACSISExtensions( const obsData * obsinfo, subSystem * subsys, unsigned int size,
		       int * status ) {
  unsigned int i;
  char type[DAT__SZTYP+1];   /* constructed type string */
  HDSLoc * temploc = NULL;
  hdsdim dim[3];
  size_t nelem;

  if (*status != SAI__OK) return;

  if (subsys->file.acsismapped) {
    *status = SAI__ERROR;
    emsRep( " ", "createACSISExtensions: ACSIS extension already mapped. Can not create", status);
    return;
  }

  /* create ACSIS extension */
  ndfXnew( subsys->file.indf, ACSISEXT, ACSISEXTTYP, 0, NULL, &(subsys->file.acsisloc), status );

  /* Need to create the following components:
     - RECEPTORS  _CHAR* array for each of the nrecep elements and their names. This fixed
       once written.
     - RECEPPOS   _DOUBLE (2 * size)   x and y positions (in tracking coordinates) for each
       receptor. This array grows in the same way as JCMTSTATE.
  */

  /* Create the receptor component and store the names */
  if (obsinfo->recep_name_buff != NULL) {
    datCctyp( obsinfo->receplen, type );
    dim[0] = obsinfo->nrecep;
    datNew( (subsys->file.acsisloc), "RECEPTORS", type, 1, dim, status );
    datFind( (subsys->file.acsisloc), "RECEPTORS", &temploc, status );
    datPut( temploc, type, 1, dim, obsinfo->recep_name_buff, status);
    datAnnul( &temploc, status );
  }

  /* Now create the positions array and map it */
  dim[0] = 2;
  dim[1] = obsinfo->nrecep;
  dim[2] = size;
  datNew( subsys->file.acsisloc, "RECEPPOS", "_DOUBLE", 3, dim, status );
  datFind( subsys->file.acsisloc, "RECEPPOS", &(subsys->file.receppos_loc), status );
  datMapD( subsys->file.receppos_loc, "WRITE", 3, dim, &(subsys->tdata.receppos), status );

  if (*status != SAI__OK) {
    subsys->file.acsismapped = 0;
    subsys->tdata.receppos = NULL;
  } else {

    /* Copy in bad values - unroll since we know there is always a 2 */
    nelem = dim[0] * dim[1] * dim[2];
    for (i=0; i < nelem; i+=2) {
      (subsys->tdata.receppos)[i] = VAL__BADD;
      (subsys->tdata.receppos)[i+1] = VAL__BADD;
    }

    subsys->file.acsismapped = 1;
  }
  
  if (*status != SAI__OK)
    emsRep(" ", "Error creating ACSIS extension", status );

  return;
}

/*
  Resize the ACSIS RECEPPOS extensions to the supplied value.
  If remap is false, the arrays will not be remapped (so call at end to resize
  before annulling locators) 
*/

static void
resizeACSISExtensions( subSystem * subsys, unsigned int newsize, 
		       int remap, int * status ) {

  hdsdim dim[3];
  size_t ndim = 3;
  int actdim;
  unsigned int j;
  unsigned int oldsize;
  unsigned int ndoubles;
  unsigned int start;

  if (*status != SAI__OK) return;

  if (subsys->file.acsismapped)
    datUnmap( subsys->file.receppos_loc, status );

  /* Get the current bounds */
  datShape( subsys->file.receppos_loc, ndim, dim, &actdim, status );

  if (*status != SAI__OK && ndim != (size_t)actdim) {
    *status = SAI__ERROR;
    emsSeti( "AD", actdim);
    emsSetu( "ND", (unsigned long) ndim );
    emsRep(" ", "Dims mismatch in ACSIS extension. ^AD != ^ND", status);
  }

  /* resize */
  oldsize = dim[ndim-1];
  dim[ndim-1] = newsize;
  datAlter( subsys->file.receppos_loc, ndim, dim, status);

  if (remap) {
    /* remap - assume this should be done after resizing all */
    datMapD( subsys->file.receppos_loc, "WRITE",
	     ndim, dim, &(subsys->tdata.receppos), status );

    if (*status == SAI__OK) {
      /* fill with bad values since may not get all receptors */
      ndoubles = dim[0] * (newsize-oldsize) * dim[1];
      start = calcOffset( dim[0], dim[1], 0, oldsize, status );
      if (*status == SAI__OK) {
	for (j=0; j < ndoubles; j++) {
	  (subsys->tdata.receppos)[start+j] = VAL__BADD;
	}
      }
    }

  }

  if (*status != SAI__OK)
    emsRep(" ", "Error resizing ACSIS extension", status );

}

/* Close down the ACSIS extension and free resources */

static void closeACSISExtensions( subSystem * subsys, int * status ) {

  if ( *status != SAI__OK ) return;

  if (subsys->curpos > 0) {
    resizeACSISExtensions( subsys, subsys->curpos, 0, status );
  }

  /* Free locators */
  datAnnul( &(subsys->file.receppos_loc), status );
  subsys->tdata.receppos = NULL;

  /* delete the receptor positions if never written */
  if (subsys->curpos == 0) {
    datErase(subsys->file.acsisloc, "RECEPPOS", status );
  }

  /* Close extension */
  datAnnul( &(subsys->file.acsisloc), status );

  subsys->file.acsismapped = 0;

  if (*status != SAI__OK)
    emsRep(" ", "Error closing ACSIS extension", status );

}

/* Write coordinate positions to ACSIS extension */
static void writeRecepPos( const obsData * obsinfo, double * posdata, unsigned int frame, 
			   const ACSISRtsState * record, int * status ) {
  unsigned int offset;

  if (*status != SAI__OK) return;

  /* Calculate offset into data array */
  offset = calcOffset( 2, obsinfo->nrecep, record->acs_feed, frame, status );

  if (posdata != NULL) {
    posdata[offset] = record->acs_feedx;
    posdata[offset+1] = record->acs_feedy;
  } else {
    *status = SAI__ERROR;
    emsRep( " ", "Attempted to write receptor positions but no data array available",
	    status );
  }

}

/* Close NDF for a particular subsystem */

static void
closeNDF( subSystem * subsys, int * status ) {

  int itemp;                /* Temp integer */
  int lbnd[NDIMS];          /* Lower bounds of NDF */
  int ubnd[NDIMS];          /* upper bounds of NDF */
  fileInfo * file;          /* File information */

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Get the file information */
  file = &(subsys->file);

  /* check that we have a file */
  if (file->indf == NDF__NOID) {
    *status = SAI__ERROR;
    emsRep("closeNDF", "attempt to close an NDF that is not actually open",
	   status);
    return;
  }

  /* Unmap */
#if SPW_DEBUG
  printf("Unmap current NDF to close the file\n");
#endif
  TIMEME("Final unmap", ndfUnmap( file->indf, "DATA", status ););

  /* Shrink file to actual size */
  ndfBound(file->indf, NDIMS, lbnd, ubnd, &itemp, status );
  if (subsys->curpos > 0) {
    ubnd[TDIM] = lbnd[TDIM] + subsys->curpos - 1;
  } else {
    ubnd[TDIM] = lbnd[TDIM];
  }

#if SPW_DEBUG
  printf("Setting final bounds. Resize to %lld spectra\n", (unsigned long long)ubnd[TDIM]);
#endif
  TIMEME( "Final set bounds", ndfSbnd(NDIMS, lbnd, ubnd, file->indf, status ););

  /* Close extensions */
  TIMEME( "Final extension resize", closeExtensions( subsys, status ););
  closeACSISExtensions( subsys, status );

  /* Close file */
  ndfAnnul( &(file->indf), status );

#if SPW_DEBUG
  printf("Wrote %d sequence steps to subsystem %d (max was %d)\n", subsys->curpos, subsys->index,
	 subsys->cursize);
#endif

  /* Force globals to be reset */
  subsys->tdata.spectra = NULL;
  subsys->cursize = 0;
  file->indf = NDF__NOID;

  if (*status != SAI__OK) {
    emsSetu("S", subsys->index);
    emsRep( " ", "Error closing Spectrum NDF file subsystem ^S", status );
  }

}


/* Resize a specific NDF */

static void
resizeNDF( const obsData * obsinfo, subSystem * subsys, unsigned int newsize, int * status ) {

  int ubnd[NDIMS];
  int lbnd[NDIMS];
  unsigned int nchans;
  unsigned int nreceps;
  unsigned int newt;
  int itemp;
  void *datapntrs[] = { NULL };/* Array of mapped pointers for ndfMap */
  size_t nbytes;
  unsigned int ncells;
  unsigned int offset;
  float * pos;
  unsigned int i;

  /* Unmap the data array */
#if SPW_DEBUG
  printf("Unmap data array in preparation for resize\n");
#endif
  TIMEME( "ndfUnmap", ndfUnmap(subsys->file.indf, "DATA", status ););

  /* Get the existing bounds */
  ndfBound(subsys->file.indf, NDIMS, lbnd, ubnd, &itemp, status );

  if (*status == SAI__OK && itemp != NDIMS) {
    *status = SAI__ERROR;
    emsSeti("N", itemp);
    emsSeti("ND", NDIMS);
    emsRep(" ", "acsSpecWriteTS: Bizarre internal error. Ndims is ^N not ^ND",
	   status);
  }
    
  nchans = ubnd[CHANDIM] - lbnd[CHANDIM] + 1;
  if (*status == SAI__OK && nchans != subsys->nchans) {
    *status = SAI__ERROR;
    emsSetu("UB", nchans);
    emsSetu("NC", subsys->nchans );
    emsRep(" ", "acsSpecWriteTS: Bizzare internal error. Nchans is ^UB not ^NC",
	   status);
  }

  nreceps = ubnd[RECDIM] - lbnd[RECDIM] + 1;
  if (*status == SAI__OK && nreceps != obsinfo->nrecep) {
    *status = SAI__ERROR;
    emsSetu("UB", nreceps);
    emsSetu("NR", obsinfo->nrecep);
    emsRep(" ", "acsSpecWriteTS: Bizzare internal error. Nreceptors is ^UB not ^NR",
	   status);
  }

  /* increment */
  ubnd[TDIM] += newsize;
  newt = ubnd[TDIM] - lbnd[TDIM] + 1;

  /* set new bounds */
#if SPW_DEBUG
  printf("Setting new bounds. Grow to %lld sequence steps (from %lld)\n", (unsigned long long)newt,
	 (unsigned long long)(newt-newsize));
#endif
  TIMEME("Set NDF bounds", ndfSbnd( NDIMS, lbnd, ubnd, subsys->file.indf, status ););

  /* map data array again */
#if SPW_DEBUG
  printf("Remap the data array\n");
#endif
  TIMEME("Remp NDF", ndfMap( subsys->file.indf, "DATA", "_REAL", "WRITE", datapntrs, &itemp, status ););
  subsys->tdata.spectra = datapntrs[0];

  /* Initialise the new memory to bad */
  ncells = ( ubnd[RECDIM] - lbnd[RECDIM] + 1 ) *
    (ubnd[CHANDIM] - lbnd[CHANDIM] + 1);
  nbytes = ncells * SIZEOF_FLOAT;
  offset = calcOffset( (ubnd[CHANDIM]-lbnd[CHANDIM]+1),
		       (ubnd[RECDIM]-lbnd[RECDIM]+1),
		       0, (newt-newsize), status); 
  if (*status == SAI__OK) {
    pos = &((subsys->tdata.spectra)[offset]);
    for (i = 0; i < newsize; i++) {
      memcpy( pos, subsys->tdata.bad, nbytes);
      pos += ncells;
    }
  }

  /* Resize the extensions */
  TIMEME("Resize extensions", resizeExtensions( subsys, newt, 1, status  ););
  TIMEME("Resize coords", resizeACSISExtensions( subsys, newt, 1, status  ););
  
  /* Update cursize */
  subsys->cursize = newt;

}

/* Calculate the offset into the 3d data array */
/* Can be used for nchan * nrecep * t data array. Use zero indexing for nrecep and tindex. */

static unsigned int
calcOffset( unsigned int nchans, unsigned int maxreceps, unsigned int nrecep, unsigned int tindex,
	    int *status ) {

  if (*status != SAI__OK) return 0;

  return (nchans * ( maxreceps * tindex + nrecep )); 

}


/* Allocate resources for spectrum */

static void
allocResources( const obsData * obsinfo, subSystem * subsys, unsigned int nseq, int *status ) {


  unsigned int seq;
  unsigned int ncells;
  size_t nbytes;
  float * pos;
  unsigned int i;


  if (*status != SAI__OK) return;

#if USE_MEMORY_CACHE
  seq = ( nseq == 0 ? subsys->maxsize : nseq );

  if (subsys->cursize != seq) {
    nbytes = seq * obsinfo->nrecep * subsys->nchans * SIZEOF_FLOAT;
    myRealloc( (void**)&(subsys->tdata.spectra), nbytes, status );

    allocHeaders(subsys, seq, status );
    allocPosData(obsinfo, subsys, seq, status );
  }

#else
  openNDF( obsinfo, subsys, subsys, nseq, status );
#endif

  /* Record the new size */
  subsys->cursize = seq;

  /* initialise all the spectra to bad in case we miss a receptor 
     in the sequence */

  pos = subsys->tdata.spectra;

  ncells = obsinfo->nrecep * subsys->nchans;
  nbytes = ncells * SIZEOF_FLOAT;
  if (*status == SAI__OK) {
    for (i = 0; i < subsys->cursize ; i++) {
      memcpy( pos, subsys->tdata.bad, nbytes );
      pos += ncells;
    }
  }

}

static void
resizeResources( const obsData *obsinfo, subSystem * subsys, unsigned int newsize, int * status ) {
  if (*status != SAI__OK) return;

#if USE_MEMORY_CACHE
  /* Currently a bug since the memory should not be realloced */
  *status = SAI__ERROR;
  emsRep(" ", "Should never be requested to realloc global buffer. Internal programming error",
	 status );
#else
  resizeNDF( obsinfo, subsys, newsize, status );
#endif

}

static void
flushResources( const obsData * obsinfo, subSystem * subsys, int * status ) {

#if USE_MEMORY_CACHE
  subSystem output;  /* Some where to store file information */

  /* initialise the output */
  memset( &output, 0, sizeof(subSystem));

  /* copy subscan number */
  output.file.subscan = subsys->file.subscan;

  /* open the file and copy in data */
  openNDF( obsinfo, subsys, &output, subsys->curpos, status );
  copyCache( obsinfo, subsys, &output, subsys->curpos, status );

  /* copy the subscan number back into the input subsystem */
  subsys->file.subscan = output.file.subscan;

  /* make sure we close the correct subsystem */
  subsys = &output;

#endif

  closeNDF( subsys, status );

  /* New position in buffer is the beginning */
  subsys->curpos = 0;

}

static void freeResources ( obsData * obsinfo, subSystem * subsys, int * status) {

  unsigned int i;

  /* Do not use status since we want to free the memory */

  if ( subsys->file.indf == NDF__NOID) {

    if ( subsys->tdata.spectra != NULL) {
      starFree( subsys->tdata.spectra);
      subsys->tdata.spectra = NULL;
    }
    for (i=0; i<NEXTENSIONS; i++) {
      if ((subsys->tdata.jcmtstate)[i] != NULL) {
	starFree( (subsys->tdata.jcmtstate)[i] );
	(subsys->tdata.jcmtstate)[i] = NULL;
      }
    }
    if (subsys->tdata.receppos != NULL) {
      starFree( subsys->tdata.receppos );
      subsys->tdata.receppos = NULL;
    }

  }

  if (subsys->tdata.bad != NULL) {
    starFree( subsys->tdata.bad );
    subsys->tdata.bad = NULL;
  }

  if (obsinfo->recep_name_buff != NULL) {
    starFree( obsinfo->recep_name_buff );
    obsinfo->recep_name_buff = NULL;
    obsinfo->receplen = 0;
  }

}

/* 
   Create space for the header information of specified size.
*/

static void
allocHeaders( subSystem * subsys, unsigned int size, int * status ) {

  int j;
  specData * tdata;

  if (*status != SAI__OK) return;

  /* Get local pointer */
  tdata = &(subsys->tdata);

  /* Make sure all are initialised */
  for (j=0; j < NEXTENSIONS; j++) {
    (tdata->jcmtstate)[j] = NULL;
  }

  /* Loop and create  */
  for (j=0; j < NEXTENSIONS; j++ ) {

    if (*status != SAI__OK) break;
    
    /* get some memory */
    /* Do not need to initialise since we always populate it */
    myRealloc( &((tdata->jcmtstate)[j]), size * hdsRecordSizes[j], status );

  }

  if (*status != SAI__OK) emsRep( " ", "allocHeaders: Unable to malloc memory for headers", status );

}

/* Allocate memory for the positional data */

static void allocPosData( const obsData * obsinfo, subSystem * subsys, unsigned int nseq, int * status ) {
  size_t nbytes;
  unsigned int ndoubles;
  unsigned int i;

  if (*status != SAI__OK) return;

  ndoubles = 2 * obsinfo->nrecep * nseq;
  nbytes = ndoubles * SIZEOF_DOUBLE;
  myRealloc( (void**)&(subsys->tdata.receppos), nbytes, status );

  /* Initialise since we can not guarantee that each receptor will get a coordinate */
  for (i=0; i<ndoubles; i++) {
    (subsys->tdata.receppos)[i] = VAL__BADD;
  }

}

/* re-allocate memory and set status - *pntr must be valid or NULL */
static void myRealloc( void ** pntr, size_t nbytes, int * status ) {
  void * tmpp;

  if (*status != SAI__OK) return;

  tmpp = starRealloc( *pntr, nbytes );
  if (tmpp == NULL) {
    *status = SAI__ERROR;
    emsSyser( "MSG", errno );
    emsSeti64( "NB", (int64_t)nbytes );
    emsRep( " ", "Unable to allocate ^NB bytes - ^MSG", status );
    starFree( *pntr );
    *pntr = NULL;
  } else {
    *pntr = tmpp;
  }

}

/* size of each type */

static size_t sizeofHDSType( const char * type, int * status ) {

  size_t retval = 0;

  switch( type[1] ) {
  case 'D':
    /* HDS sizes are known in advance */
    retval = 8;
    break;
  case 'I':
  case 'R':
    retval = 4;
    break;
  case 'C':
    /* offset past the _CHAR* */
    retval = strtol( &(type[6]) , NULL, 10);
    if (retval == 0) {
      *status = SAI__ERROR;
      emsSyser( "MESSAGE", errno );
      emsSetc( "TYP", type );
      emsRep(" ", "Unable to determine length of string '^TYP' - ^MESSAGE", status );
    }
    break;
  default:
    *status = SAI__ERROR;
    emsSetc( "TYP", type );
    emsRep( " ", "Error determining size of supplied type '^TYP'", status );
  }

  return retval;

}

/* Copy from malloced buffers to mapped buffers */

 static void copyCache( const obsData * obsinfo, const subSystem * input, subSystem * output, 
		       unsigned int nseq, int * status) {
  unsigned int i;

  if (*status != SAI__OK) return;

  /* check compatibility */
  if (input->nchans != output->nchans) {
    *status = SAI__ERROR;
    emsSetu( "IN", input->nchans);
    emsSetu( "OUT", output->nchans);
    emsRep( "copyCache", "Number of channels in input subsystem (^IN) differs from"
	    " number in output subsystem (^OUT)", status);
    return;
  }

  if (output->maxsize < nseq) {
    *status = SAI__ERROR;
    emsSetu( "NSEQ", nseq);
    emsSetu( "MAX", output->maxsize);
    emsRep("copyCache2", "Number of sequence steps to copy (^NSEQ) exceeds available"
	   " space in output (^MAX)", status);
    return;
  }

  /* Copy the main data array */
  memcpy( output->tdata.spectra, input->tdata.spectra, obsinfo->nrecep * input->nchans * nseq * SIZEOF_FLOAT );
  output->curpos = nseq;

  /* sequence data */
  for (i=0; i<NEXTENSIONS; i++) {
    memcpy( (output->tdata.jcmtstate)[i], (input->tdata.jcmtstate)[i], nseq * hdsRecordSizes[i] );
  }

  /* receptor positions */
  memcpy( output->tdata.receppos, input->tdata.receppos, obsinfo->nrecep * 2 * nseq * SIZEOF_DOUBLE );

}



/********************************** Debug functions ********************/

#if SPW_DEBUG
/* simply subtract two timeval structs and return the answer */
/* Does tp2 - tp1 */
static double duration ( struct timeval * tp1, struct timeval * tp2 ) {
  double diff = 0.0;
  diff = (tp2->tv_sec - tp1->tv_sec) +
    (tp2->tv_usec - tp1->tv_usec ) / 1E6;

  return diff;
}
#endif
