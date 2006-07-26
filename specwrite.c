
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
#include <unistd.h>
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

/* Application name for history writing */
#define APPNAME "ACSIS-DA (V" PACKAGE_VERSION ")"

/* Enable memory cache */
#define USE_MEMORY_CACHE 1

/* pi/180:  degrees to radians */
#define DD2R 0.017453292519943295769236907684886127134428718885417

/* Debug prints
   0 - disabled
   1 - standard debug 
   2 - verbose debug
 */
#define SPW_DEBUG_LEVEL 0

/* Largest file name allowed (including path) */
#define MAXFILE 1024

/* Define the number of extensions we support */
#define NEXTENSIONS 45

/* Maximum number of subsystems we can handle 
   We know that ACSIS can have at most 4 spectral subsystems
   and they will not change during a single observation. */

#define MAXSUBSYS 4
static const unsigned int maxsubsys = MAXSUBSYS;

/* String size for obsData */
#define SPW__SZFSTAT 10

/* Global state variables */

/* This struct gives an overview of the observation state itself */
typedef struct obsData {
  /* number of receptors in this observation */
  unsigned int nrecep;
  size_t receplen;       /* longest receptor name */
  char * recep_name_buff;  /* buffer for receptor names */
  char focal_station[SPW__SZFSTAT+1];  /* Focal station of this instrument */
  float * fplanex;        /* X and Y coordinates of receptors in focal plane arcsec */
  float * fplaney;
  /* observation number, ut date and number of subsystems */
  unsigned int obsnum;
  unsigned int yyyymmdd;
  unsigned int nsubsys; /* Actual number of subsystems in use */
  /* relevant data directories */
  char datadir[MAXFILE+1];   /* current data directory */
  char rootdir[MAXFILE+1];   /* root directory */
} obsData;

/* actual data (either in memory or mapped from disk) */
typedef struct specData {
  float  * spectra;  /* Array of data to receive spectra ( nchans x nrecep x nseq) */
  double * receppos; /* Receptor positions (2 x  nrecep x nseq) */
  float  * tsys;     /* System temp (nrecep * nseq) */
  void   * jcmtstate[NEXTENSIONS];  /* Pointers to JCMTSTATE information */
  float  * bad;      /* array of bad values to easily initialise new time slice */
  unsigned char * count;    /* array of size (nrecep x nseq) containing 1 if spectrum written
                               else 0. Can be used to indicate when all spectra received
                               or if overwriting a spectrum. */
} specData;

/* This struct contains file information (ndf identifiers, hds locators) */
typedef struct fileInfo {
  unsigned int subscan;  /* current subscan number */
  int indf;   /* NDF identifier for this file */
  HDSLoc * extloc;  /* JCMSTATE extension locator */
  HDSLoc * acsisloc; /* ACSIS extension locator */
  HDSLoc * extlocators[NEXTENSIONS]; /* Locators to each JCMTSTATE component */
  HDSLoc * receppos_loc;   /* Locators to each .MORE.ACSIS.RECEPPOS */  
  HDSLoc * tsys_loc;       /* Locator to .MORE.ACSIS.TSYS */
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
  int          inseq;   /* We are in a sequence */
  unsigned int seqlen;  /* Expected length of this sequence */
  int          alloced; /* called allocResources no this struct? */
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
static char * getFileName( const char * dir, unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status );
static char * getOkFileName( const char * dir, unsigned int yyyymmdd,
			   unsigned int obsnum, int * status );
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

static void mapThisExtension( HDSLoc * loc, size_t ndim, unsigned int oldtsize, unsigned int newtsize,
		       const char type[], void ** mapped, int * status );

static void resizeThisExtension ( HDSLoc * loc, size_t ndim, unsigned int newtsize,
				  int ismapped, unsigned int * oldtsize, int * status );

static void writeRecord( void * basepntr[], unsigned int tindex,
			 const ACSISRtsState * record,
			 int * status );
static void writeRecepPos( const obsData * obsinfo, double * posdata, unsigned int tindex,
			   const ACSISRtsState * record, int * status );
static void writeTSys( const obsData * obsinfo, float * data, unsigned int frame, 
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
static void allocTsysData( const obsData * obsinfo, subSystem * subsys, unsigned int nseq, 
			   int * status );
static void
flushResources( const obsData * obsinfo, subSystem * subsys, int * status );

static void copyCache( const obsData * obsinfo, const subSystem * input,
		       subSystem * output, unsigned int nseq, int * status);

static size_t sizeofHDSType( const char * type, int * status );

static void myRealloc( void **pntr, size_t nbytes, int * status );

static int hasAllSpectra( const obsData * obsinfo, const subSystem * subsys,
			  int * status );
static int hasSeqSpectra( const obsData * obsinfo, const subSystem * subsys,
			  unsigned int tindex, int * status );

void writeFlagFile (const obsData * obsinfo, const subSystem subsystems[],
		    int * status);

void writeWCSandFITS (const obsData * obsinfo, const subSystem subsystems[],
		      const AstFitsChan * fits[], int * status);

AstFrameSet *specWcs( const AstFrameSet *fs, int ntime, const double times[], int * status );

static void checkNoFileExists( const char * file, int * status );

#if SPW_DEBUG_LEVEL > 0
static double duration ( struct timeval * tp1, struct timeval * tp2 );
#endif

/* Stolen code */
static int kpgPtfts( int indf, const AstFitsChan * fchan, int * status );


/* Function to put quotes around a symbol so that we can do
   CPP string concatenation */
#define myxstr(s) mystr(s)
#define mystr(s) #s
#define CHARTYP(s) "_CHAR*" myxstr(s)

/* Threshold for reporting timing anomalise */
#define LONGTIME 2.0

/* Macro to time an event */
#if SPW_DEBUG_LEVEL > 0
#define TIMEME(label,func) { struct timeval tp1; struct timeval tp2; double tpdiff; \
    gettimeofday( &tp1, NULL ); \
    func;			\
    gettimeofday( &tp2, NULL ); \
    tpdiff = duration( &tp1, &tp2 ); \
    if (tpdiff > LONGTIME) printf( ">>>>>>>>>>>>>>>>>" label " took %.3f seconds <<<<<<<<<<<<<<<\n", tpdiff); \
  }
#else
#define TIMEME(label, func)  func
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
#define RTS_NUM           0
#define RTS_END           1
#define RTS_TASKS         2
#define SMU_X             3
#define SMU_Y             4
#define SMU_Z             5
#define SMU_CHOP_PHASE    6
#define SMU_JIG_INDEX     7
#define SMU_AZ_JIG_X      8
#define SMU_AZ_JIG_Y      9
#define SMU_AZ_CHOP_X    10
#define SMU_AZ_CHOP_Y    11
#define SMU_TR_JIG_X     12
#define SMU_TR_JIG_Y     13
#define SMU_TR_CHOP_X    14
#define SMU_TR_CHOP_Y    15
#define TCS_AIRMASS      16
#define TCS_AZ_ANG       17
#define TCS_AZ_AC1       18
#define TCS_AZ_AC2       19
#define TCS_AZ_DC1       20
#define TCS_AZ_DC2       21
#define TCS_AZ_BC1       22
#define TCS_AZ_BC2       23
#define TCS_BEAM         24
#define TCS_INDEX        25
#define TCS_SOURCE       26
#define TCS_TR_SYS       27
#define TCS_TR_ANG       28
#define TCS_TR_AC1       29
#define TCS_TR_AC2       30
#define TCS_TR_DC1       31
#define TCS_TR_DC2       32
#define TCS_TR_BC1       33
#define TCS_TR_BC2       34
#define JOS_DRCONTROL    35
#define ENVIRO_REL_HUM   36
#define ENVIRO_PRESSURE  37
#define ENVIRO_AIR_TEMP  38
#define POL_ANG          39
#define ACS_SOURCE_RO    40
#define ACS_NO_PREV_REF  41
#define ACS_NO_NEXT_REF  42
#define ACS_NO_ONS       43
#define ACS_EXPOSURE     44

/* Definitions of HDS types associated with ACSISRtsStates struct. All these
   will be created in the file. */
static const char * hdsRecordNames[NEXTENSIONS][2] = 
  {
   { "_INTEGER", "RTS_NUM" },
   { "_DOUBLE", "RTS_END" },
   { CHARTYP(SIZEOF_RTS_TASKS), "RTS_TASKS" },
   { "_DOUBLE", "SMU_X" },
   { "_DOUBLE", "SMU_Y" },
   { "_DOUBLE", "SMU_Z" },
   { CHARTYP(SIZEOF_SMU_CHOP_PHASE), "SMU_CHOP_PHASE" },
   { "_INTEGER", "SMU_JIG_INDEX" },
   { "_DOUBLE", "SMU_AZ_JIG_X" },
   { "_DOUBLE", "SMU_AZ_JIG_Y" },
   { "_DOUBLE", "SMU_AZ_CHOP_X" },
   { "_DOUBLE", "SMU_AZ_CHOP_Y" },
   { "_DOUBLE", "SMU_TR_JIG_X" },
   { "_DOUBLE", "SMU_TR_JIG_Y" },
   { "_DOUBLE", "SMU_TR_CHOP_X" },
   { "_DOUBLE", "SMU_TR_CHOP_Y" },
   { "_DOUBLE", "TCS_AIRMASS" },
   { "_DOUBLE", "TCS_AZ_ANG" },
   { "_DOUBLE", "TCS_AZ_AC1" },
   { "_DOUBLE", "TCS_AZ_AC2" },
   { "_DOUBLE", "TCS_AZ_DC1" },
   { "_DOUBLE", "TCS_AZ_DC2" },
   { "_DOUBLE", "TCS_AZ_BC1" },
   { "_DOUBLE", "TCS_AZ_BC2" },
   { CHARTYP(SIZEOF_TCS_BEAM), "TCS_BEAM" },
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
   { "_INTEGER", "JOS_DRCONTROL" },
   { "_REAL", "ENVIRO_REL_HUM" },
   { "_REAL", "ENVIRO_PRESSURE" },
   { "_REAL", "ENVIRO_AIR_TEMP" },
   { "_DOUBLE", "POL_ANG" },
   { CHARTYP(SIZEOF_ACS_SOURCE_RO), "ACS_SOURCE_RO" },
   { "_INTEGER", "ACS_NO_PREV_REF" },
   { "_INTEGER", "ACS_NO_NEXT_REF" },
   { "_INTEGER", "ACS_NO_ONS" },
   { "_REAL", "ACS_EXPOSURE" }
  };

/* Somewhere to store the precomputed sizes of each HDS element */
static size_t hdsRecordSizes[NEXTENSIONS];

/* Extension support */

/* Name of STATE and ACSIS extensions - some could almost be shared with SCUBA2... */
#define STATEEXT   "JCMTSTATE"
#define STATEEXTTYPE "RTS_ARR"
#define ACSISEXT   "ACSIS"
#define ACSISEXTTYP "ACSIS_COMP"

/* Mandatory FITS headers written by spec writer itself */
/* Name of the FITS header containing the subscan information. */

#define FITS_NSUBSCAN "NSUBSCAN"

/* Name of the FITS header indicating when the last subscan is found */

#define FITS_OBSEND "OBSEND"

/* Name of fits header for temperature scale */

#define FITS_TEMPSCAL "TEMPSCAL"

/*********************** NDF "cube" FILE *************************************/

/* Number of sequences to increment file size by if it runs out of space */

/* May want this parameter to be settable as the number of spectra that
   we expect in a given time period */
#define MAXRECEP   16
#define MAXRATE    20
#define PRESIZETIME 10
#define NGROW  (MAXRATE * PRESIZETIME)

/* Number of bytes we should accumulate before opening a new file */
#define MAXBYTES ( 512 * 1024 * 1024 )

/* if we have unfeasibly small spectra and 1 receptor we could get
   a very large requirement for memory for all the extensions
   associated with MAXBYTES above. Therefore limit the number of 
   sequence steps we can actual get.
   Peak rate is 10 MB/s which would be maxseq of 
          MAXBYTES / (4 * 8192 * MAXRECEP) = 1024 sequences
   whereas min rate gives
          MAXBYTES / (4 *    1 * 1       ) = 134 million sequences
   or more feasibly
          MAXBYTES / (4 * 1024 * 1       ) = 130,000 sequences
 */
#define MAXSEQ ( MAXBYTES / ( 4 * 1024 * 1 ) )

/* maximum number of sequence steps we can get out of sequence */
#define MAXSEQERR  5


/*
*+
*  Name:
*     acsSpecOpenTS

*  Purpose:
*     Prepare the file writing system for a new observation

*  Invocation:
*     acsSpecOpenTS( const char * dir, unsigned int yyyymmdd, 
*                    unsigned int obsnum, unsigned int nrecep,
*                    unsigned int nsubsys, const char *recepnames[],
*                    const char * focal_station,
*                    const float fplanex[], const float fplaney[],
*                    int * status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to prepare the data writing system.
*     It must be called before calling acsSpecWriteTS. The subscan
*     directory will be created.

*  Arguments:
*     dir = const char * (Given)
*        Directory to write the file.
*     yyyymmdd = unsigned int (Given)
*        UT date in YYYYMMDD format. Used to construct the file name.
*     obsnum = unsigned int (Given)
*        Current observation number. Used to construct the file name.
*        There is no check to make sure that the observation number is new.
*     nrecep = unsigned int (Given)
*        Number of receptors participating in this observation.
*     nsubsys = unsigned int (Given)
*        Number of subsystems present in this observation.
*     recepnames[] = const char*[] (Given)
*        Names of each receptor in feed order.
*     focal_station = const char* (Given)
*        Focal station for the instrument. [DIRECT, NASMYTH_L, NASMYTH_R]
*     fplanex[] = const float[] (Given)
*        X offsets of each receptor in the focal plane (arcsec).
*     fplaney[] = const float[] (Given)
*        Y offsets of each receptor in the focal plane (arcsec).
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.
*     05-APR-2006 (TIMJ):
*        Use structured globals
*     21-APR-2006 (TIMJ):
*        Defer resource allocation until the first spectrum arrives.
*     26-JUL-2005 (TIMJ):
*        Add FPLANE arguments.

*  Notes:
*     - Currently only one observation can be active at any time
*       (some global variables are used internally for state)
*     - It is an error for this routine to be called if an observation
*       is in progress. Please call acsSpecCloseTS() before calling
*       acsSpecOpenTS().
*     - The files created by this routine conform to the ICD (OCS/ICD/022)

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
	       const char *recepnames[],
	       const char * focal_station,
	       const float fplanex[], const float fplaney[],
	       int * status ) {

  char * cpos = NULL;          /* offset into string */
  unsigned int i;              /* Loop counter */
  char * sdir = NULL;          /* subscan directory */
  size_t len;                  /* temp length */
  size_t receplen;             /* Length of longest receptor name */
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
  if (focal_station != NULL) {
    strncpy( OBSINFO.focal_station, focal_station, SPW__SZFSTAT );
    (OBSINFO.focal_station)[SPW__SZFSTAT] = '\0';
  } else {
    (OBSINFO.focal_station)[0] = '\0';
  }
  OBSINFO.fplanex = NULL;
  OBSINFO.fplaney = NULL;

  /* Create the directory to receive each subscan */
  sdir = createSubScanDir( dir, yyyymmdd, obsnum, status );

  /* Store the resulting directory and root dir */
  if (*status == SAI__OK) {
    strncpy(OBSINFO.datadir, sdir, MAXFILE);
    (OBSINFO.datadir)[MAXFILE] = '\0';
    strncpy(OBSINFO.rootdir, dir, MAXFILE);
    (OBSINFO.rootdir)[MAXFILE] = '\0';


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

    /* Allocate memory for the focal plane offsets */
    if (*status == SAI__OK && fplanex != NULL && fplaney != NULL) {
      OBSINFO.fplanex = starMalloc( nrecep * sizeof(*fplanex) );
      OBSINFO.fplaney = starMalloc( nrecep * sizeof(*fplaney) );
      memcpy( OBSINFO.fplanex, fplanex, nrecep * sizeof(*fplanex) );
      memcpy( OBSINFO.fplaney, fplaney, nrecep * sizeof(*fplaney) );
    }

    /* Store dynamic receptor information */
    OBSINFO.receplen = receplen;
    OBSINFO.recep_name_buff = recep_name_buff;

    /* Need an NDF/cache per subsystem. Make sure everything is initialised here.
       We allocate the resources the first time that acsSpecWriteTS is called. */
    for (i = 0; i < nsubsys; i++) {

      /* Select the subsystem to modify */
      subsys = &(SUBSYS[i]);

      /* zero it out */
      memset( subsys, 0, sizeof(*subsys));

      /* Some intialisation */
      subsys->index = i;

    }

    /* indicate that an observation is now ready */
    INPROGRESS = 1;

  }

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
*     result = acsSpecWriteTS( unsigned int subsys, unsigned int nchans, 
*                     const float spectrum[], const ACSISRtsState * record,
*                     int *status);

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to write a spectrum to an HDS container
*     that had previously been opened using acsSpecOpen.

*  Arguments:
*     subsys = unsigned int (Given)
*        Subsystem used for this spectrum (start counting at 0).
*     nchans = unsigned int (Given)
*        Number of channels in the spectrum. The number of channels
*        must be identical for all spectra supplied for the identical
*        subsystem.
*     spectrum = float[nchan] (Given)
*        Spectrum itself.
*     record = const ACSISRtsState * (Given)
*        Header information associated with this spectrum.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  Returned Value:
*     acsSpecWriteTS = int
*        Returns 1 if the spectrum was accepted as a science
*        spectrum. Returns -1 if the spectrum was processed as
*        a calibration. Value is undefined if bad status is returned.

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.
*     05-APR-2006 (TIMJ):
*        Use structured globals
*     21-APR-2006 (TIMJ):
*        Dynamically allocate memory based on first spectrum arriving.
*     09-MAY-2006 (TIMJ):
*        Keep track of all spectra arriving so we can determine
*        when all spectra for a sequence have arrived.
*     11-MAY-2006 (TIMJ):
*        Remove Freq argument.

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

int
acsSpecWriteTS( unsigned int subsysnum, unsigned int nchans, const float spectrum[], 
		const ACSISRtsState* record,
	        int * status ) {

  float * data; /* local copy of mapped pointer to spectrum */
  unsigned int offset;         /* offset into data array */
  unsigned int coff;           /* offset into count array */
  int seqinc = 0;              /* did we increment sequence number? */
  unsigned int ngrow;          /* number of time slices to grow file */
  unsigned int reqnum;         /* number of time slices indicates by RTS sequence */
  unsigned int tindex;         /* Position in sequence array for this sequence */
  unsigned int *rtsseqs;       /* Pointer to array of sequence numbers */
  unsigned int max_behind;     /* How many sequence steps to look behind */
  unsigned int maxseq;         /* calculated maximum number of sequence steps allowed */
  unsigned int nperseq;        /* Number of elements per sequence */
  unsigned int inpos;          /* curpos on entry */
  unsigned int last_seqnum = 0; /* Sequence number encountered for previous spectrum in this slot */

  int found = 0;               /* Did we find the sequence? */
  unsigned int i;              /* loop counter */
  double * posdata;            /* pointer to position data */
  float  * tsysdata;           /* pointer to Tsys data */
  void ** recdata;
  unsigned int startind;
  subSystem * subsys;

  ACSISRtsState  state;        /* local editable copy of state information */

  if (*status != SAI__OK) return 0;

  /* make sure that the subsys number is in range */
  if ( subsysnum >= maxsubsys ) {
    *status = SAI__ERROR;
    emsSetu("IN", subsysnum);
    emsSeti("MAX", maxsubsys-1);
    emsRep(" ","acsSpecWriteTS: Supplied subsystem number (^IN) exceeds max allowed (^MAX)", status);
    return 0;
  }
  if (subsysnum >= OBSINFO.nsubsys) {
    *status = SAI__ERROR;
    emsSetu( "IN", subsysnum );
    emsSetu( "MAX", OBSINFO.nsubsys -1 );
    emsRep( " ", "acsSpecWriteTS: Supplied subsystem number (^IN) exceeds number supplied to "
	    "acsSpecOpenTS (^MAX)", status);
    return 0;
  }

  /* Check to see if we've already been called */
  if (!INPROGRESS) {
    *status = SAI__ERROR;
    emsRep("HDS_SPEC_WRITETS_ERR1",
	   "acsSpecWriteTS called, yet an observation has not been initialised", status);
    return 0;
  }

  /* Check feed range */
  if ( record->acs_feed >= OBSINFO.nrecep ) {
    *status = SAI__ERROR;
    emsSetu( "NR", OBSINFO.nrecep - 1 );
    emsSetu( "FEED", record->acs_feed );
    emsRep( " ", "acsSpecWriteTS called, yet the feed number (^FEED) exceeds the expected number (^NR)",
	    status );
  }

  /* check that we are SPECTRUM_RESULT. In the future this will trigger the
     use of a CALDATA extension for calibrations but that is not yet
     implemented. */
  if ( strncmp( record->acs_source_ro, "SPECTRUM_RESULT", SIZEOF_ACS_SOURCE_RO )
       != 0) {
    return -1;
  }

  /* first need to take a local copy for editing */
  memcpy( &state, record, sizeof(ACSISRtsState) );

  /* Some elements are compatibility elements for SCUBA-2 and should be
     cleared here so that ACSIS does not have to worry about them */
  strcpy( state.smu_chop_phase, " " );
  strcpy( state.tcs_beam, " " );
  state.smu_az_chop_x = VAL__BADD;
  state.smu_az_chop_y = VAL__BADD;
  state.smu_tr_chop_x = VAL__BADD;
  state.smu_tr_chop_y = VAL__BADD;

  /* also convert the feed coordinates from degrees to radians */
  state.acs_feedx *= DD2R;
  state.acs_feedy *= DD2R;

  /* Get local copy of subsystem from global */
  subsys = &(SUBSYS[subsysnum]);

  /* store the input curpos */
  inpos = subsys->curpos;

  /* The reference number of channels depends on whether this spectrum is meant
     to go in the SPECTRUM_RESULT or CALDATA part of the file (since CALDATA
     can have more channels). */

  /* if this is the first time through, we need to allocate resources based on the number
     of channels that we have been given. */

  if (subsys->nchans == 0) {

    /* Number of data values per sequence */
    nperseq = OBSINFO.nrecep * nchans;

    /* Store the number of channels for this subsystem */
    subsys->nchans = nchans;

    /* Calculate the number of sequence steps that we are allowed to grow
       before opening new file. */
    maxseq = MAXBYTES / ( nperseq * SIZEOF_FLOAT );
    subsys->maxsize = ( maxseq > MAXSEQ ? MAXSEQ : maxseq );

#if SPW_DEBUG_LEVEL > 1
    printf("Calculated maximum number of sequence steps per file: %u\n", subsys->maxsize);
#endif

    /* Allocate a cache of bad values to simplify initialisation */
    subsys->tdata.bad = starMalloc( nperseq * SIZEOF_FLOAT );
    if (subsys->tdata.bad == NULL) {
      if (*status == SAI__OK) {
	*status = SAI__ERROR;
	emsRep(" ","Unable to allocate memory for bad value cache", status );
      }
    } else {
      for (i=0; i<nperseq; i++) {
	(subsys->tdata.bad)[i] = VAL__BADR;
      }
    }

    /* Allocate the count array */
    subsys->tdata.count = starMalloc( OBSINFO.nrecep * subsys->maxsize * sizeof(unsigned char) );
    if (subsys->tdata.count == NULL) {
      if (*status == SAI__OK) {
	*status = SAI__ERROR;
	emsRep(" ","Unable to allocate memory for count array", status );
      }
    } else {
      memset( subsys->tdata.count, 0, OBSINFO.nrecep * subsys->maxsize * sizeof(unsigned char) );
    }

  } else if ( subsys->nchans != nchans ) {
    *status = SAI__ERROR;
    emsSetu( "IN", nchans );
    emsSetu( "REF", subsys->nchans );
    emsSetu( "SS", subsysnum );
    emsRep( " ", "acsSpecWriteTS: Number of channels in subsystem ^SS has changed "
	    "from ^REF to ^IN", status );
    return 0;
  }

  /* Allocate resources for this subsystem if not currently allocated */
  if (!subsys->alloced) {
#if SPW_DEBUG_LEVEL > 1
    printf("+++++++ Need to allocate resources on entry to WriteTS\n");
#endif
    allocResources( &OBSINFO, subsys, subsys->maxsize, status );
  }

  /* first thing to do is determine whether this sequence number is new or old */

  if (subsys->curpos == 0) {
    /* have not written anything yet so the correct place to write this
       sequence is at index 0 in the data array (assuming cursize[] is
       large enough) */
    seqinc = 1;
    tindex = 0;
    /* printf("First spectrum for this file for this subsystem\n"); */

  } else {
    /* if the supplied value is the most recent value then we do not 
       need to search */

    if (state.rts_num == subsys->curseq) {

      tindex = subsys->curpos - 1;
      /* printf("Reusing sequence %u at index %u\n", subsys->curseq, tindex); */

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
	if (state.rts_num == rtsseqs[startind-i]) {
	  /* found the sequence */
	  tindex = startind - i;
	  found = 1;
	  break;
	}
      }

      if (found == 0) {
#if SPW_DEBUG_LEVEL > 1
	printf("Did not find sequence %u\n", state.rts_num);
#endif
	/* did not find this sequence number so it is a new one */
	tindex = subsys->curpos; /* curpos will be incremented */
	seqinc = 1;
      } else {
	/* going back in time so this may not be efficient */
#if SPW_DEBUG_LEVEL > 1
	printf("Sequence %u matches index %u. Previous seq num=%u\n", 
	       state.rts_num, tindex, subsys->curseq);
#endif
	subsys->curseq = state.rts_num;

      }
    }
  }
  
  /* if the sequence number has incremented we need to increase the t-axis counter */

  if (seqinc) {

#if SPW_DEBUG_LEVEL > 1
    printf("Incrementing sequence number to sequence %u\n", state.rts_num);
#endif

    /* store the new value */
    subsys->curseq = state.rts_num;

    /* increment the counters value */
    (subsys->curpos)++;

    /* See if we need to grow */
    /* We can grow either because we have suddenly realised we don't fit *or* because
       we have been told how many sequence steps to expect - calculate the required
       number to extend. (but we know at least 1) */
    reqnum = 1;

    if ( state.rts_endnum > state.rts_num ) {
      reqnum = state.rts_endnum - state.rts_num + 1;
    }

    /* Calculate the length of this sequence if it is started by this
       sequence step. */
    if (!subsys->inseq) {
#if SPW_DEBUG_LEVEL > 1
      printf("+++++++++++++++++++++SEQLEN set to %u\n",reqnum);
#endif
      subsys->seqlen = reqnum;
      subsys->inseq = 1; /* we are now in a sequence */
    }


    /* if the required number exceeds the maximum allowed size then 
       we need to reduce the reqnum to a more manageable level and
       hope that we can flush in the middle of a sequence. We could simply
       set it to one more than the current size (we just need to prevent
       the file from being flushed when empty) but if resizing is required,
       we may as well give a good hint rather than resize each time around.
     */
    if ( reqnum > (subsys->maxsize - subsys->curpos) ) {
#if SPW_DEBUG_LEVEL > 1
      printf("Length of sequence = %u Maxspace = %u Currently=%u\n",
	      reqnum, subsys->maxsize, subsys->curpos);
#endif
      reqnum = subsys->maxsize - subsys->curpos;
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
#if SPW_DEBUG_LEVEL > 1
	printf("--------flush after unexpected grow\n");
#endif
	flushResources( &OBSINFO, subsys, status );

#if USE_MEMORY_CACHE
	/* always want to make sure we allocate the max amount of memory */
	ngrow = subsys->maxsize;
#endif

#if SPW_DEBUG_LEVEL > 1
	printf("------- alloc after unexpected grow\n");
#endif
	allocResources( &OBSINFO, subsys, ngrow, status );

	/* indicate that we are starting at the beginning with the next spectrum */
	tindex = 0;

	/* current position must be 1 now since we still have to write the
	   spectrum, so force it to that.
	    - allocResources does not set this, but flushResources resets to 0 */
	subsys->curpos = 1;

      } else {
#if SPW_DEBUG_LEVEL > 0
	printf("Cursize: %u Curpos: %u ngrow: %u maxsize: %u ; Need to resize.\n",
	       subsys->cursize, subsys->curpos, ngrow, subsys->maxsize); 
#endif
	/* Resize the NDF */
	resizeResources( &OBSINFO, subsys, ngrow, status );
      }
    }
  }

  /* copy in the data */
  if (*status == SAI__OK) {

#if SPW_DEBUG_LEVEL > 0
    if (seqinc) {
      printf(">>> About to write first spectrum at position %u for sequence %u\n", tindex,
	     subsys->curseq);
    }
#endif
    /* Calculate offset into array - number of spectra into the array times number of
       channels per spectrum. */
    offset = calcOffset( subsys->nchans, OBSINFO.nrecep,
			 state.acs_feed, tindex, status );

    /* Get local copies of pointers */
    data = (subsys->tdata.spectra);
    recdata = (subsys->tdata.jcmtstate);
    posdata = (subsys->tdata.receppos);
    tsysdata = (subsys->tdata.tsys);

    /* We would like to know the sequence number that was used previously for this
       time slot */
    last_seqnum = 0;
    if (!seqinc) last_seqnum = (((unsigned int *)recdata[RTS_NUM])[tindex]);

    /* Calculate offset into count array: a  */
    coff = calcOffset( 1, OBSINFO.nrecep, state.acs_feed, tindex, status );

    /* check to make sure this slot is free */
    if ( (subsys->tdata.count)[coff] != 0) {
      if (*status == SAI__OK) {
	*status = SAI__ERROR;
	emsSetu("CURSEQ", subsys->curseq );
	emsSetu("PREVSEQ", last_seqnum);
	emsSetu("FEED", state.acs_feed );
	if ( (subsys->tdata.count)[coff] == 1 ) {
	  emsRep(" ", "acsSpecWriteTS: Error. Overwriting a slot that already contains a spectrum"
		 " (current sequence number = ^CURSEQ, previous sequence was ^PREVSEQ,"
		 " feed number = ^FEED)", status );
	} else {
	  emsSetu("NWRITE", (unsigned int)(subsys->tdata.count)[coff]);
	  emsRep( " ","acsSpecWriteTS: Error. Bizarrely have already written ^NWRITE"
		  " spectra to this slot (current sequence number = ^CURSEQ, feed = ^FEED)",
		  status);
	}
      }
    }

    /* sanity check the sequence numbers in this spectrum with that already stored
       if we have not just started a new sequence */
    if ( !seqinc && subsys->curseq != last_seqnum) {
      if (*status == SAI__OK) {
	*status = SAI__ERROR;
	emsSetu( "LAST", last_seqnum );
	emsSetu( "CUR", subsys->curseq);
	emsRep(" ","Last time this slot was used it had sequence number ^LAST"
	       " but now it has value ^CUR.", status);
      }
    }

#if SPW_DEBUG_LEVEL > 1
    printf(">><<>> Writing spectrum from feed %u to tindex %u curpos %u offset %u\n",
	   state.acs_feed, tindex, subsys->curpos, offset);
#endif

    if (*status == SAI__OK)
      memcpy( &(data[offset]), spectrum, subsys->nchans*SIZEOF_FLOAT );

    /* Store record data and receptor positions. Base record only updates each
       sequence step but receptor positions and tsys should be written for all records. */
    if (seqinc) writeRecord( recdata, tindex, &state, status );
    writeRecepPos( &OBSINFO, posdata, tindex, &state, status );
    writeTSys( &OBSINFO, tsysdata, tindex, &state, status );

    /* increment the count for this location */
    if (*status == SAI__OK) (subsys->tdata.count)[coff]++;

    /* see whether this sequence is complete */
    if ( hasSeqSpectra( &OBSINFO, subsys, tindex, status ) ) {

#if SPW_DEBUG_LEVEL > 0
      printf("<<< Sequence %u (tindex=%u) completed with this spectrum for feed %u\n",
	     subsys->curseq, tindex, state.acs_feed);
#endif
      /* We *could* flush to disk at this point if we do not think the next sequence
	 will fit. This may alleviate network contention if we end up dumping the
	 file at the start of a sequence. This will only trigger if the current
	 sequence is complete. We may need to check that the previous few are
	 also complete if we get them in random order.
      */
      if ( subsys->curseq == state.rts_endnum ) {

	/* Sequence complete so no longer in a sequence */
	subsys->inseq = 0;

	/* are we now complete for all sequences so far? */
	if (hasAllSpectra( &OBSINFO, subsys, status) ) {

	  /* will we have space for the next sequence? */
	  if ( (subsys->maxsize - subsys->curpos) < subsys->seqlen ) {
	    /* write what we have pre-emptively even if the next sequence
	       turns out to be shorter */
#if SPW_DEBUG_LEVEL > 0
	    printf("!!!!!!! COMPLETE SEQUENCE & FULL BUFFER <<<<<<<<\n");
#endif
	    flushResources( &OBSINFO, subsys, status );

	    /* Do not call allocResources() since we might not get
	       another spectrum */

	  }
	}
      }
    }

    /* reset insequence flag if this was the last sequence step in it*/
    if (subsys->curseq == state.rts_endnum) subsys->inseq = 0;

#if SPW_DEBUG_LEVEL > 1
    if (seqinc && *status == SAI__OK) {
      printf(">>> Wrote first spectrum "
#  if USE_MEMORY_CACHE
	     "to memory cache "
#  else
	     "to disk "
#  endif  /* USE_MEMORY_CACHE */
	     "at this position %u (curpos = %u)\n", tindex, subsys->curpos );
    }
#endif  /* SPW_DEBUG_LEVEL */

    /* Sanity check - curpos must either be 1 or one bigger than
       earlier if we added a new sequence. It can not be 0 if a spectrum
       was supplied and we added a new sequence.
    */
    if (*status == SAI__OK && seqinc) {
      if (subsys->curpos != 1 && subsys->curpos != (inpos + 1) ) {
	*status = SAI__ERROR;
	emsSetu( "IN", inpos);
	emsSetu( "OUT", subsys->curpos);
	emsRep(" ", "acsSpecWriteTS: Needed to write new sequence but seemingly overwote:"
	       " (input position was ^IN, output was ^OUT)", status);
      }
    }

  } /* status ok for writing a spectrum */

  return 1;
}

/*
*+
*  Name:
*     acsSpecCloseTS

*  Purpose:
*     Write FITS header and close HDS file.

*  Invocation:
*     acsSpecCloseTS( const AstFitsChan * fits[], int incArchiveBounds, int *status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to close the file after all spectra have been
*     written. The FITS header is written.

*  Arguments:
*     fits[] = const AstFitsChan * (Given)
*        Array of FITS headers. One per subsystem.
*     incArchiveBounds = int (Given)
*        If true, bounds FITS keywords will be calculated.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.
*     20-APR-2006 (TIMJ):
*        Use structured globals.
*     09-MAY-2006 (TIMJ):
*        Attempt to run even if status is bad on entry.
*     11-MAY-2006 (TIMJ):
*        Add flag for adding archive bounds to FITS header.

*  Notes:
*     - Must have previously called acsSpecOpenTS.
*     - File is resized to the actual number of spectra written.
*     - If status is bad on entry this routine will attempt to execute
*       so that spectra will be written rather than deleted from internal
*       memory.

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
acsSpecCloseTS( const AstFitsChan * fits[], int incArchiveBounds, int * status ) {

  unsigned int i;           /* Loop counter */
  int found = 0;            /* Found an open NDF? */
  subSystem * subsys;       /* specific subsystem */
  int lstat = SAI__OK;      /* Local status */

  /* Ideally should not automatically abort on entry since we may be required
     to save the spectra that we have already been given. The trick is in proceeding
     if we can report that spectra have been written (either because curpos is non-zero
     or because the subscan number is non-zero). If status was bad on entry we should
     not add to the error message stack. If it was good on entry we should add to the
     stack. This requires conditional use of emsMark and emsRlse and using an internal
     error status value.
  */

  /* Check to see if we've already been called */
  if (!INPROGRESS) {
    if (*status == SAI__OK) {
      *status = SAI__ERROR;
      emsRep("HDS_SPEC_CLOSETS_ERR1",
	     "acsSpecCloseTS called, yet an observation has not been initialised",
	     status);
    }
    return;
  }

  /* Start new error context if status is bad */
  if (*status != SAI__OK) emsMark();

  /* Loop over each subsystem and flush cached spectra to disk */
  found = 0;
  for (i = 0; i < OBSINFO.nsubsys; i++) {
    /* Get local copy of subsystem from global */
    subsys = &(SUBSYS[i]);
    if ( subsys->file.indf != NDF__NOID || subsys->tdata.spectra != NULL) {
      found = 1;
#if SPW_DEBUG_LEVEL > 1
      printf("-------- Final close flush\n");
#endif
      flushResources( &OBSINFO, subsys, &lstat);
    }
    freeResources( &OBSINFO, subsys, &lstat );
  }

  /* report error if not found any open NDFs */
  if (lstat == SAI__OK && !found) {
    lstat = SAI__ERROR;
    emsRep(" ", "acsSpecCloseTS: Failed to find open NDF components", &lstat );
  }

  /* Now need to open all the files that we have opened previously and adjust
     any FITS headers.
  */

#if SPW_DEBUG_LEVEL > 1
  astShow( fits[0] );
#endif

  writeWCSandFITS( &OBSINFO, SUBSYS, fits, &lstat );

  /* Write the flag file contents */

  writeFlagFile( &OBSINFO, SUBSYS, &lstat);

  /* Force globals to be reset */
  for (i = 0; i < maxsubsys; i++) {
    memset( &(SUBSYS[i]), 0, sizeof(SUBSYS[i]) );
  }
  memset( &OBSINFO, 0, sizeof(OBSINFO) );

  /* reset progress */
  INPROGRESS = 0;

  if (lstat != SAI__OK) {
    emsRep( " ", "Error closing Spectrum file", &lstat );
  }

  /* release the mark if status on entry was bad so that 
     we do not contaminate the error stack. */
  if (*status != SAI__OK) {
    if (lstat != SAI__OK) emsAnnul( &lstat );
    emsRlse();
  } else {
    /* make sure that the exit status is returned */
    *status = lstat;
  }
 
}

/*********************** HELPER FUNCTIONS (PRIVATE) *****************************/

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

/* Form the OK file name
   - returns a pointer to static memory.
 */

static char * getOkFileName( const char * dir, unsigned int yyyymmdd,
			     unsigned int obsnum, int * status ) {

  static char filename[MAXFILE]; /* buffer for filename - will be returned */
  int flen;                        /* Length of string */
  char * root;                   /* Root filename */

  if (*status != SAI__OK) return NULL;
  root = getFileRoot( yyyymmdd, 0, obsnum, 0, status );
  if (*status != SAI__OK) return NULL;

  /* Form the file name - assume posix filesystem */
  flen = snprintf(filename, MAXFILE, "%s/.%s.ok", dir, root );

  if (flen >= MAXFILE) {
    *status = SAI__ERROR;
    emsSeti("SZ", MAXFILE );
    emsSetu("N", obsnum );
    emsSetu("UT", yyyymmdd );
    emsRep("HDS_SPEC_OPEN_ERR1",
	   "Error forming Ok filename. Exceeded buffer size of ^SZ chars for scan ^N on UT ^UT", status );
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
   - if "dir" is NULL or zero length  no directory is prefixed
   - returns a pointer to static memory.
 */

static char * getDirName( const char * dir, unsigned int yyyymmdd,
			  unsigned int obsnum, int * status ) {

  static char dirname[MAXFILE]; /* buffer for dirname - will be returned */
  int flen;                        /* Length of string */

  if (*status != SAI__OK) return NULL;

  /* check for null directory  and make sure we do not have a separator
     if null */
  if ( dir != NULL && strlen(dir) > 0) {
    /* Form the file name - assume posix filesystem */
    flen = snprintf(dirname, MAXFILE, "%s/%05u", dir, obsnum );
  } else {
    flen = snprintf(dirname, MAXFILE, "%05u", obsnum );
  }

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

  checkNoFileExists( dir, status );

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
  char * history[1] = { "Initial writing of raw data." };

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

#if SPW_DEBUG_LEVEL > 0
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

#if SPW_DEBUG_LEVEL > 0
  printf("Opening NDF file '%s' to default size of %u sequence steps\n", ndfname, ngrow);
#endif

  /* Sanity check - want to make sure that we do not open a file that
     is already there. There could be a race condition between the
     check and the rename() but that should be impossible in normal DA
     usage so we ignore the possibility.
  */
  checkNoFileExists( ndfname, status );

  /* create the NDF */
  ndfPlace( NULL, ndfname, &place, status );
  ndfNew( "_REAL", NDIMS, lbnd, ubnd, &place, &(file->file.indf), status );

  /* Update the cursize[] array and the nchans array */
  file->cursize = ubnd[TDIM] - lbnd[TDIM] + 1;

  /* History component */
  ndfHcre( file->file.indf, status );
  ndfHput("NORMAL",APPNAME, 1, 1, history,
	  0, 0, 0, file->file.indf, status );

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

  /* Real Time Sequencer */
  ((int *)basepntr[RTS_NUM])[frame] = record->rts_num;
  ((double *)basepntr[RTS_END])[frame] = record->rts_end;

  cnfExprt( record->rts_tasks,
	   (char *)basepntr[RTS_TASKS]+
	    SIZEOF_RTS_TASKS*frame,
	    SIZEOF_RTS_TASKS );

  /* Secondary Mirror Unit */
  ((double *)basepntr[SMU_X])[frame] = record->smu_x;
  ((double *)basepntr[SMU_Y])[frame] = record->smu_y;
  ((double *)basepntr[SMU_Z])[frame] = record->smu_z;

  cnfExprt( record->smu_chop_phase,
	    (char *)basepntr[SMU_CHOP_PHASE]+
	    SIZEOF_SMU_CHOP_PHASE*frame,
	    SIZEOF_SMU_CHOP_PHASE);

  ((int *)basepntr[SMU_JIG_INDEX])[frame] = record->smu_jig_index;
  ((double *)basepntr[SMU_AZ_JIG_X])[frame] = record->smu_az_jig_x;
  ((double *)basepntr[SMU_AZ_JIG_Y])[frame] = record->smu_az_jig_y;
  ((double *)basepntr[SMU_AZ_CHOP_X])[frame] = record->smu_az_chop_x;
  ((double *)basepntr[SMU_AZ_CHOP_Y])[frame] = record->smu_az_chop_y;
  ((double *)basepntr[SMU_TR_JIG_X])[frame] = record->smu_tr_jig_x;
  ((double *)basepntr[SMU_TR_JIG_Y])[frame] = record->smu_tr_jig_y;
  ((double *)basepntr[SMU_TR_CHOP_X])[frame] = record->smu_tr_chop_x;
  ((double *)basepntr[SMU_TR_CHOP_Y])[frame] = record->smu_tr_chop_y;

  /* Telescope Control System */
  ((double *)basepntr[TCS_AIRMASS])[frame] = record->tcs_airmass;
  ((double *)basepntr[TCS_AZ_ANG])[frame] = record->tcs_az_ang;
  ((double *)basepntr[TCS_AZ_AC1])[frame] = record->tcs_az_ac1;
  ((double *)basepntr[TCS_AZ_AC2])[frame] = record->tcs_az_ac2;
  ((double *)basepntr[TCS_AZ_DC1])[frame] = record->tcs_az_dc1;
  ((double *)basepntr[TCS_AZ_DC2])[frame] = record->tcs_az_dc2;
  ((double *)basepntr[TCS_AZ_BC1])[frame] = record->tcs_az_bc1;
  ((double *)basepntr[TCS_AZ_BC2])[frame] = record->tcs_az_bc2;

  cnfExprt( record->tcs_beam,
	    (char *)basepntr[TCS_BEAM]+
	    SIZEOF_TCS_BEAM*frame,
	    SIZEOF_TCS_BEAM);

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

  /* JOS control */
  ((int *)basepntr[JOS_DRCONTROL])[frame] = record->jos_drcontrol;

  /* ENVIRO task */
  ((float *)basepntr[ENVIRO_AIR_TEMP])[frame] = record->enviro_air_temp;
  ((float *)basepntr[ENVIRO_PRESSURE])[frame] = record->enviro_pressure;
  ((float *)basepntr[ENVIRO_REL_HUM])[frame] = record->enviro_rel_hum;

  /* POLarimeter aka ROVER */
  ((double *)basepntr[POL_ANG])[frame] = record->pol_ang;

  /* ACSIS internal */
  cnfExprt( record->acs_source_ro,
	    (char*)basepntr[ACS_SOURCE_RO]+ SIZEOF_ACS_SOURCE_RO*frame,
	    SIZEOF_ACS_SOURCE_RO );

  ((int *)basepntr[ACS_NO_PREV_REF])[frame] = record->acs_no_prev_ref;
  ((int *)basepntr[ACS_NO_NEXT_REF])[frame] = record->acs_no_next_ref;
  ((int *)basepntr[ACS_NO_ONS])[frame] = record->acs_no_ons;
  ((float *)basepntr[ACS_EXPOSURE])[frame] = record->acs_exposure;

}

/* Create the .MORE.ACSIS extensions (that are not JCMT state structure members) */


static void
createACSISExtensions( const obsData * obsinfo, subSystem * subsys, unsigned int size,
		       int * status ) {
  char type[DAT__SZTYP+1];   /* constructed type string */
  HDSLoc * temploc = NULL;
  hdsdim dim[3];
  void * tpntr = NULL;

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
     - RECEPPOS   _DOUBLE (2 * nrecep * size)   x and y positions 
       (in tracking coordinates) for each
       receptor. This array grows in the same way as JCMTSTATE.
     - TSYS       _REAL (nrecep * size) Grows as JCMTSTATE grows.
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

  /* Create the FOCAL_STATION component and store it */
  if (obsinfo->nrecep > 0 && strlen(obsinfo->focal_station) > 0 ) {
    datNew0C( (subsys->file.acsisloc), "FOCAL_STATION", strlen(obsinfo->focal_station),
	      status);
    datFind( (subsys->file.acsisloc), "FOCAL_STATION", &temploc, status );
    datPut0C( temploc, obsinfo->focal_station, status);
    datAnnul( &temploc, status );
  }

  /* Focal plane positions */
  if ( obsinfo->fplanex != NULL ) {
    datNew1R( (subsys->file.acsisloc), "FPLANEX", obsinfo->nrecep, status );
    datFind( (subsys->file.acsisloc), "FPLANEX", &temploc, status );
    datPut1R( temploc, obsinfo->nrecep, obsinfo->fplanex, status);
    datAnnul( &temploc, status );
  }
  if ( obsinfo->fplaney != NULL ) {
    datNew1R( (subsys->file.acsisloc), "FPLANEY", obsinfo->nrecep, status );
    datFind( (subsys->file.acsisloc), "FPLANEY", &temploc, status );
    datPut1R( temploc, obsinfo->nrecep, obsinfo->fplaney, status);
    datAnnul( &temploc, status );
  }

  /* Now create the positions array and map it */
  dim[0] = 2;
  dim[1] = obsinfo->nrecep;
  dim[2] = size;
  datNew( subsys->file.acsisloc, "RECEPPOS", "_DOUBLE", 3, dim, status );
  datFind( subsys->file.acsisloc, "RECEPPOS", &(subsys->file.receppos_loc), status );
  mapThisExtension( subsys->file.receppos_loc, 3, 0, size, "_DOUBLE",
		    &tpntr, status);
  subsys->tdata.receppos = tpntr;
  tpntr = NULL;

  /*  datMapD( subsys->file.receppos_loc, "WRITE", 3, dim, &(subsys->tdata.receppos), status );
   */

  /* Now the TSYS array and map it */
  dim[0] = obsinfo->nrecep;
  dim[1] = size;
  datNew( subsys->file.acsisloc, "TSYS", "_REAL", 2, dim, status );
  datFind( subsys->file.acsisloc, "TSYS", &(subsys->file.tsys_loc), status);
  mapThisExtension( subsys->file.tsys_loc, 2, 0, size, "_REAL",
		    &tpntr, status );
  subsys->tdata.tsys = tpntr;

  /*  datMapR( subsys->file.tsys_loc, "WRITE", 2, dim, &(subsys->tdata.tsys),
	   status );
  */
  if (*status != SAI__OK) {
    subsys->file.acsismapped = 0;
    subsys->tdata.receppos = NULL;
    subsys->tdata.tsys = NULL;
  } else {
    subsys->file.acsismapped = 1;
  }
  
  if (*status != SAI__OK)
    emsRep(" ", "Error creating ACSIS extension", status );

  return;
}

/*
  Resize the ACSIS RECEPPOS and TSYS extensions to the supplied value.
  If remap is false, the arrays will not be remapped (so call at end to resize
  before annulling locators) 
*/

static void
resizeACSISExtensions( subSystem * subsys, unsigned int newsize, 
		       int remap, int * status ) {

  unsigned int old_rpos_size;
  unsigned int old_tsys_size;
  void * tpntr = NULL;

  if (*status != SAI__OK) return;

  /* RECEPPOS */
  resizeThisExtension( subsys->file.receppos_loc, 3, newsize,
		       subsys->file.acsismapped, &old_rpos_size, status );

  /* TSYS */
  resizeThisExtension( subsys->file.tsys_loc, 2, newsize,
		       subsys->file.acsismapped, &old_tsys_size, status );


  /* update mapped status */
  subsys->file.acsismapped = 0;

  /* Now remap these extensions */
  if (remap && *status == SAI__OK) {
    /* in principal old_rpos_size should equal old_tsys_size */
    mapThisExtension( subsys->file.receppos_loc, 3, old_rpos_size, newsize, "_DOUBLE",
		      &tpntr, status );
    subsys->tdata.receppos = tpntr;
    tpntr = NULL;
    mapThisExtension( subsys->file.tsys_loc, 2, old_tsys_size, newsize, "_REAL",
		      &tpntr, status );
    subsys->tdata.tsys = tpntr;
    subsys->file.acsismapped = ( *status == SAI__OK ? 1 : 0 );
  }

  if (*status != SAI__OK)
    emsRep(" ", "Error resizing ACSIS extension", status );

}

/* resize the extension - unmapping first if necessary */

static void resizeThisExtension ( HDSLoc * loc, size_t ndim, unsigned int newtsize,
				  int ismapped, unsigned int * oldtsize, int * status ) {

  hdsdim dim[DAT__MXDIM];
  int actdim;

  *oldtsize = 0;

  if (*status != SAI__OK) return;

  if (ismapped) {
    datUnmap( loc, status );
  }

  /* Get the current bounds */
  datShape( loc, ndim, dim, &actdim, status );

  if (*status == SAI__OK && ndim != (size_t)actdim) {
    *status = SAI__ERROR;
    emsSeti( "AD", actdim);
    emsSetu( "ND", (unsigned int) ndim );
    emsRep(" ", "Dims mismatch in ACSIS extension during resizing. ^AD != ^ND", status);
  }

  /* resize */
  if (*status == SAI__OK) {
    *oldtsize = dim[ndim-1];
    dim[ndim-1] = newtsize;
    datAlter( loc, ndim, dim, status);
  }
}

/* Map the extension but also fills with bad values */

static void mapThisExtension( HDSLoc * loc, size_t ndim, unsigned int oldtsize, unsigned int newtsize,
		       const char type[], void ** mapped, int * status ) {

  hdsdim dim[DAT__MXDIM];
  int actdim;
  unsigned int j;
  unsigned int nelems = 0;
  unsigned int start;
  unsigned int secdim = 0;  /* second dimension */

  if ( *status != SAI__OK) return;

  /* Get the current bounds */
  datShape( loc, ndim, dim, &actdim, status );

  if (*status == SAI__OK && ndim != (size_t)actdim) {
    *status = SAI__ERROR;
    emsSeti( "AD", actdim);
    emsSetu( "ND", (unsigned int) ndim );
    emsRep(" ", "Dims mismatch in ACSIS extension during mapping. ^AD != ^ND", status);
  }

  /* map */
  datMap( loc, type, "WRITE", ndim, dim, mapped, status );

  /* Calculate how many elements we need to fill in */
  if ( ndim == 3) {
    nelems = dim[0] * (newtsize-oldtsize) * dim[1];
    secdim = dim[1];
  } else if (ndim == 2) {
    nelems = dim[0] * (newtsize-oldtsize);
    secdim = 1; /* hack so that calcOffset works in 2D */
  } else {
    if ( *status == SAI__OK ) {
      *status = SAI__ERROR;
      emsSetu( "ND",  ndim);
      emsRep( " ", "Extension mapping routine can only handle 2 and 3 dimensions not ^ND - internal programming error. ",
	      status );
    }
  }

  /* only fill if we have new elements */
  if (*status == SAI__OK && nelems > 0) {
    /* fill with bad values since may not get all receptors */
    start = calcOffset( dim[0], secdim, 0, oldtsize, status );

    if (*status == SAI__OK) {

      /* need to switch on type */
      if ( strcmp( type, "_DOUBLE") == 0 ) {
	for (j=0; j < nelems; j++) {
	  ((double *)*mapped)[start+j] = VAL__BADD;
	}
      } else if ( strcmp( type, "_REAL" ) == 0 ) {
	for (j=0; j < nelems; j++) {
	  ((float *)*mapped)[start+j] = VAL__BADR;
	}
      } else {
	*status = SAI__ERROR;
	emsSetc( "TYP", type );
	emsRep(" ", "Unrecognized data type: ^TYP - internal programming error",
	       status);
      }
    }

  }


}

/* Close down the ACSIS extension and free resources */

static void closeACSISExtensions( subSystem * subsys, int * status ) {

  if ( *status != SAI__OK ) return;

  if (subsys->curpos > 0) {
    resizeACSISExtensions( subsys, subsys->curpos, 0, status );
  }

  /* Free locators */
  datUnmap( subsys->file.receppos_loc, status );
  datAnnul( &(subsys->file.receppos_loc), status );
  subsys->tdata.receppos = NULL;
  datUnmap( subsys->file.tsys_loc, status );
  datAnnul( &(subsys->file.tsys_loc), status );
  subsys->tdata.tsys = NULL;

  /* delete the receptor positions if never written */
  if (subsys->curpos == 0) {
    datErase(subsys->file.acsisloc, "RECEPPOS", status );
    datErase(subsys->file.acsisloc, "TSYS", status );
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

/* Write Tsys to ACSIS extension */
static void writeTSys( const obsData * obsinfo, float * data, unsigned int frame, 
		       const ACSISRtsState * record, int * status ) {
  unsigned int offset;

  if (*status != SAI__OK) return;

  /* Calculate offset into data array */
  offset = calcOffset( 1, obsinfo->nrecep, record->acs_feed, frame, status );

  if (data != NULL) {
    data[offset] = record->acs_tsys;
  } else {
    *status = SAI__ERROR;
    emsRep( " ", "Attempted to write Tsys information but no data array available",
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
#if SPW_DEBUG_LEVEL > 0
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

#if SPW_DEBUG_LEVEL > 0
  printf("Setting final bounds. Resize to %lld time steps\n", (unsigned long long)ubnd[TDIM]);
#endif
  TIMEME( "Final set bounds", ndfSbnd(NDIMS, lbnd, ubnd, file->indf, status ););

  /* Close extensions */
  TIMEME( "Final extension resize", closeExtensions( subsys, status ););
  closeACSISExtensions( subsys, status );

  /* Close file */
  ndfAnnul( &(file->indf), status );

#if SPW_DEBUG_LEVEL > 0
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
#if SPW_DEBUG_LEVEL > 0
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
#if SPW_DEBUG_LEVEL > 0
  printf("Setting new bounds. Grow to %lld sequence steps (from %lld)\n", (unsigned long long)newt,
	 (unsigned long long)(newt-newsize));
#endif
  TIMEME("Set NDF bounds", ndfSbnd( NDIMS, lbnd, ubnd, subsys->file.indf, status ););

  /* map data array again */
#if SPW_DEBUG_LEVEL > 0
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


/* Allocate resources for spectrum and initialise */

static void
allocResources( const obsData * obsinfo, subSystem * subsys, unsigned int nseq, int *status ) {


  unsigned int seq;
  unsigned int ncells;
  size_t nbytes;
  float * pos;
  unsigned int i;


  if (*status != SAI__OK) return;

  seq = ( nseq == 0 ? subsys->maxsize : nseq );

#if USE_MEMORY_CACHE

  if (subsys->cursize != seq) {
    nbytes = seq * obsinfo->nrecep * subsys->nchans * SIZEOF_FLOAT;
    myRealloc( (void**)&(subsys->tdata.spectra), nbytes, status );

    allocHeaders(subsys, seq, status );
    allocPosData(obsinfo, subsys, seq, status );
    allocTsysData(obsinfo, subsys, seq, status );
  }

#else
  openNDF( obsinfo, subsys, subsys, nseq, status );
#endif

  /* count array is always in memory */
  nbytes = obsinfo->nrecep * seq * sizeof( unsigned char );
  myRealloc( (void**)&(subsys->tdata.count), nbytes, status );
  if (*status == SAI__OK) memset( subsys->tdata.count, 0, nbytes );

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

  /* indicate that we have been allocated */
  subsys->alloced = 1;

}

static void
resizeResources( const obsData *obsinfo, subSystem * subsys, unsigned int newsize, int * status ) {
  if (*status != SAI__OK) return;

#if USE_MEMORY_CACHE
  /* Currently a bug since the memory should not be realloced */
  *status = SAI__ERROR;
  emsSetu("SZ", newsize);
  emsSetu("N", subsys->index);
  emsSetu("OBS", obsinfo->obsnum);
  emsRep(" ", "Should never be requested to realloc global buffer by ^SZ sequences"
	 " (subsys ^N, obs ^OBS)."
	 " Internal programming error",
	 status );
#else
  resizeNDF( obsinfo, subsys, newsize, status );
#endif

}

/* Write current state to file or close existing file. Does not clear
   struct contents. */

static void
flushResources( const obsData * obsinfo, subSystem * subsys, int * status ) {

  subSystem * toclose;  /* pointer to subsystem struct we are actually closing */
#if SPW_DEBUG_LEVEL > 0
  double percent = 0.0;
#endif

  if (*status != SAI__OK) return;

#if USE_MEMORY_CACHE
  subSystem output;  /* Some where to store file information */

#if SPW_DEBUG_LEVEL > 0
  if (subsys->cursize > 0) {
    percent = 100.0 * (double)subsys->curpos / (double)subsys->cursize;
  }
  printf("Flushing with memory cache = %u/%u (%.1f%% capacity)\n", subsys->curpos,
	 subsys->cursize, percent);
#endif /* SPW_DEBUG_LEVEL */

  /* if we have no spectra we have nothing to write so do not want
     to open the NDF - should we set status? */
  if (subsys->curpos == 0) return;

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
  toclose = &output;
#else
  toclose = &subsys;

#endif  /* USE_MEMORY_CACHE */

  closeNDF( toclose, status );

  /* New position in buffer is the beginning */
  subsys->curpos = 0;
  subsys->curseq = 0;

  /* indicate that we need to clear this resource */
  subsys->alloced = 0;

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
    if (subsys->tdata.tsys != NULL) {
      starFree( subsys->tdata.tsys );
      subsys->tdata.tsys = NULL;
    }

  }

  if (subsys->tdata.bad != NULL) {
    starFree( subsys->tdata.bad );
    subsys->tdata.bad = NULL;
  }

  if (subsys->tdata.count != NULL) {
    starFree( subsys->tdata.count );
    subsys->tdata.count = NULL;
  }

  if (obsinfo->recep_name_buff != NULL) {
    starFree( obsinfo->recep_name_buff );
    obsinfo->recep_name_buff = NULL;
    obsinfo->receplen = 0;
  }

  if (obsinfo->fplanex != NULL) {
    starFree( obsinfo->fplanex );
    obsinfo->fplanex = NULL;
  }

  if (obsinfo->fplaney != NULL) {
    starFree( obsinfo->fplaney );
    obsinfo->fplaney = NULL;
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

/* Allocate memory for the Tsys data */

static void allocTsysData( const obsData * obsinfo, subSystem * subsys, unsigned int nseq, int * status ) {
  size_t nbytes;
  unsigned int nelems;
  unsigned int i;

  if (*status != SAI__OK) return;

  nelems = 2 * obsinfo->nrecep * nseq;
  nbytes = nelems * SIZEOF_FLOAT;
  myRealloc( (void**)&(subsys->tdata.tsys), nbytes, status );

  /* Initialise since we can not guarantee that each receptor will get a value */
  for (i=0; i<nelems; i++) {
    (subsys->tdata.tsys)[i] = VAL__BADR;
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

  /* TSys */
  memcpy( output->tdata.tsys, input->tdata.tsys, obsinfo->nrecep * nseq * SIZEOF_FLOAT );

}

 /* See whether the "count" array has all slots filled for this time slice */

static int hasSeqSpectra( const obsData * obsinfo, const subSystem * subsys,
			  unsigned int tindex, int * status ) {

  unsigned int offset; /* offset into count data array */
  unsigned int i;    /* loop counter */
  int missing = 0; /* true if we found a missing spectrum */

  if (*status != SAI__OK) return 0;

  offset = calcOffset( 1, obsinfo->nrecep, 0, tindex, status );

  if (*status == SAI__OK) {
    for (i=0; i < obsinfo->nrecep; i++) {
      if ( (subsys->tdata.count)[offset+i] == 0 ) {
	missing = 1;
	break;
      }
    }
  }

  /* return true is missing is false */
  return ( missing ? 0 : 1 );
}

/* See whether all spectra have been received for all sequences
   received so far. */

static int hasAllSpectra( const obsData * obsinfo, const subSystem * subsys,
			  int * status ) {
  unsigned int last; /* end position into count data array */
  unsigned int i;    /* loop counter */
  int missing = 0; /* true if we found a missing spectrum */

  if (*status != SAI__OK) return 0;

  /* get last valid offset in array */
  last = calcOffset( 1, obsinfo->nrecep, obsinfo->nrecep - 1, 
		     subsys->curpos - 1, status );
  if (*status == SAI__OK) {
    for (i=0; i < last; i++) {
      if ( (subsys->tdata.count)[i] == 0 ) {
	missing = 1;
	break;
      }
    }
  }

  /* return true is missing is false */
  return ( missing ? 0 : 1 );
}
			  
/* Write the OK flag file and content */

void writeFlagFile (const obsData * obsinfo, const subSystem subsystems[],
		    int * status) {

  int fd;                   /* file descriptor of temp ok file */
  FILE * fstream = NULL;    /* Stream associated with temp ok file */
  char * fname;             /* file name of a given subscan */
  int flen;                 /* length of return value from snprintf */
  unsigned int i;           /* Loop counter */
  unsigned int j;           /* Loop counter */
  char * ldir = NULL;       /* subscan directory relative to rootdir */
  char tmpok[MAXFILE];      /* temporary file name */
  char *okfile;             /* name of ok file */
  int sysstat;              /* system status return */
  const subSystem * subsys;       /* specific subsystem */

  if (*status != SAI__OK) return;

  /* We must first write to a temporary file so that we can rename it when it
     is complete. */

  flen = snprintf(tmpok, MAXFILE,
#if HAVE_MKSTEMPS
		  "%s/tempXXXXXX.ok",
#elif HAVE_MKSTEMP
		  "%s/tempXXXXXX",
#endif
		  obsinfo->rootdir );
  if (flen >= MAXFILE && status == SAI__OK) {
    *status = SAI__ERROR;
    emsRep(" ","Error forming temporary 'ok' filename. Exceeded buffer",
	   status );
  }

  /* get a temporary file */
  if (*status == SAI__OK) {
#if HAVE_MKSTEMPS
    fd = mkstemps( tmpok, 3 ); /* .ok is 3 characters */
#elif HAVE_MKSTEMP
    fd = mkstemp( tmpok );
#else
    UNABLE TO CREATE A TEMPORARY FILE
#endif
    if (fd == -1) {
      *status = SAI__ERROR;
      emsSyser( "ERRNO", errno );
      emsSetc( "DIR", obsinfo->rootdir );
      emsSetc( "TMP", tmpok );
      emsRep(" ","acsSpecCloseTS: Failed to open temporary 'ok' file ^TMP in dir ^DIR:"
	     " ^ERRNO", status );
    }
  }

  /* Force permissions on the temp file so that the pipeline can read it */
  if (*status == SAI__OK) {
    sysstat = fchmod(fd, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    if (sysstat == -1) {
	*status = SAI__ERROR;
	emsSyser( "ERRNO", errno );
	emsRep(" ","acsSpecCloseTS: Error setting permissions on 'ok' file: ^ERRNO",
	       status );
    }
  }

  /* use stream I/O so open this on a stream */
  if (*status == SAI__OK) {
    fstream = fdopen( fd, "w+" );
    if (!fstream) {
      *status = SAI__ERROR;
      emsSyser( "ERRNO", errno );
      emsRep(" ","acsSpecCloseTS: Failed to open stream on 'ok' file descriptor:"
	     " ^ERRNO", status);
    }
  }

  /* for each subsystem, for each subscan, write a FITS header and entry in ok
     file */
  if (*status == SAI__OK) {

    for (i = 0; i < obsinfo->nsubsys; i++) {
      /* local subsystem */
      subsys = &(subsystems[i]);

      /* subscans start counting at 1 */
      for (j = 1; j <= subsys->file.subscan ; j++ ) {
	/* contents of okay file are relative to data dir but should
	   not include datadir itself */
	ldir = getDirName( NULL, obsinfo->yyyymmdd, obsinfo->obsnum, status );
	fname = getFileName( ldir, obsinfo->yyyymmdd, i,
			     obsinfo->obsnum, j, status );
	fprintf( fstream, "%s\n", fname );
      }

    }
  }

  /* rename temp file to okfile and close temp file */
  okfile = getOkFileName( obsinfo->rootdir, obsinfo->yyyymmdd,
			  obsinfo->obsnum, status );

  /* Sanity check - want to make sure that we do not open a file that
     is already there. There could be a race condition between the
     check and the rename() but that should be impossible in normal DA
     usage so we ignore the possibility.
  */
  checkNoFileExists( okfile, status );

  if (*status == SAI__OK) {
    sysstat = rename(tmpok, okfile);
    if (sysstat == -1) {
	*status = SAI__ERROR;
	emsSyser( "ERRNO", errno );
	emsSetc( "F", tmpok);
	emsSetc( "T", okfile );
	emsRep(" ","acsSpecCloseTS: Error renaming okay file from '^F'"
	       " to ^T: ^ERRNO", status);
    }
  }
  /* close temporary file descriptor and stream (even if bad status) */
  if (fstream) fclose( fstream ); /* also closes fd */

  return;
}

void writeWCSandFITS (const obsData * obsinfo, const subSystem subsystems[],
		      const AstFitsChan * fits[], int * status) {

  char * fname;             /* file name of a given subscan */
  unsigned int i;           /* Loop counter */
  int          indf;        /* NDF identifier */
  unsigned int j;           /* Loop counter */
  const subSystem * subsys = NULL; /* Current subsystem */
  int          place;       /* unused NDF placeholder */
  AstFitsChan * lfits = NULL; /* Local copy of FITS header */
  int *        oldstat;     /* Internal AST status on entry */
  AstFrameSet *  wcs = NULL; /* World coordinates system from FITS header */
  HDSLoc * xloc = NULL;    /* Extension locator */
  size_t tsize;            /* number of sequence steps in file */
  HDSLoc * tloc = NULL;    /* Locator to time data */
  double * tdata = NULL;   /* Pointer to mapped MJD time data */
  int      tempscal = 0;   /* temperature scale? */
  void * tpntr = NULL;     /* temporary generic pointer */
  AstFrameSet * specwcs = NULL; /* framset for timeseries cube */
  char * stemp = NULL;     /* temporary pointer to string */
  int sysstat;             /* status from system call */
  char * history[1] = { "Finalise headers and make ICD compliant." };

  /* headers to be retained */
  char * retainfits[] = {
    "DATE-OBS", "OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z", "SSYSOBS", NULL
  };

  if (*status != SAI__OK) return;

  /* do nothing if we have null pointer */
  if (fits == NULL) return;

  /* AST routines so register status */
  oldstat = astWatch( status );

  /* Start an AST context so we do not need to annul AST pointers explicitly. */
   astBegin;

  /* loop over each sub system and open each subscan */

  /* repeat code in flag file function rather than have a separate
     routine that populates an array of strings (or even store those
     filenames in the subsys struct itself). This is just laziness */

  for (i = 0; i < obsinfo->nsubsys; i++) {
    /* local subsystem */
    subsys = &(subsystems[i]);

    /* subscans start counting at 1 */
    for (j = 1; j <= subsys->file.subscan ; j++ ) {
      /* need full path of file */
      fname = getFileName( obsinfo->datadir, obsinfo->yyyymmdd, i,
			   obsinfo->obsnum, j, status );

      /* open the NDF */
#if SPW_DEBUG_LEVEL > 0
      printf("Writing file FITS header %s\n", fname );
#endif
      ndfOpen( NULL, fname, "UPDATE", "OLD", &indf, &place, status );

      /* manipulate FITS header here...First take a copy. */
      lfits = astCopy( fits[i] );

      /* need to add a SUBSCAN number to the header */
      astClear( lfits, "Card" );
      astFindFits( lfits, FITS_NSUBSCAN, NULL, 0 );
      astSetFitsI( lfits, FITS_NSUBSCAN, (int)j, "Sub-scan number", 1);

      /* need to add a OBSEND number to the header. True if
	 this is the last file */
      astClear( lfits, "Card" );
      astFindFits( lfits, FITS_OBSEND, NULL, 0 );
      astSetFitsL( lfits, FITS_OBSEND,
		   ( j == subsys->file.subscan ? 1 : 0 ),
		   "True if file is last in current observation", 1);

      /* Need to look for the BUNIT header */
      astClear( lfits, "Card");
      tempscal = 0;
      if ( astGetFitsS( lfits, "BUNIT", &stemp ) ) {
	ndfCput( stemp, indf, "UNITS", status );
	if (strncmp("K", stemp,1) == 0) tempscal = 1;
	astDelFits( lfits );
      } else {
	/* not calibrated */
	ndfCput( "uncalibrated", indf, "UNITS", status );
      }

      /* attach a data label */
      if (tempscal) {
	/* Note the use of AST control codes for subscript/superscript */
	ndfCput( "T%s60+%v30+A%^50+%<20+*%+   corrected antenna temperature", indf, "LABEL", status );
	/* and look for a TEMPSCAL fits header  - which should be undef by default and
	   so can stay undef if we are uncalibrated */
	astClear( lfits, "Card" );
	if ( astFindFits( lfits, FITS_TEMPSCAL, NULL, 0 ) ) {
	  astSetFitsS(lfits, FITS_TEMPSCAL, "TA*", "Temperature scale in use", 1);
	}

      } else {
	ndfCput( "Power", indf, "LABEL", status );
      }
      
      /* Remove the END card */
      astClear( lfits, "Card" );
      if (astFindFits( lfits, "END", NULL, 0 ) ) {
	astDelFits( lfits );
      }

      /* Bounds associated with this file */

      /* Mark some headers to be retained after stripping */
      i = 0;
      while ( retainfits[i] != NULL ) {
	/* clearing each time is not very efficient but I don't want
	   to burn in the ordering */
	astClear( lfits, "Card" );
	astFindFits( lfits, retainfits[i], NULL, 0 );
	astRetainFits( lfits );
	i++;
      }

      /* write astrometry */
      /* Rewind the FitsChan */
      astClear( lfits, "Card" );

      /* extract astrometry from the FITS header */
      wcs = astRead( lfits );

      /* need the time information */
      ndfXloc( indf, STATEEXT, "READ", &xloc, status );
      datFind( xloc, "RTS_END", &tloc, status );
      
      datMapV( tloc, "_DOUBLE", "READ", &tpntr, &tsize, status );
      tdata = tpntr;

      /* calculate the frameset */
      specwcs = specWcs( wcs, (int)tsize, tdata, status);

      /* clean up */
      datUnmap( tloc, status );
      datAnnul( &tloc, status );
      datAnnul( &xloc, status );

      /* Write WCS */
      ndfPtwcs( specwcs, indf, status );

      /* write FITS header */
      kpgPtfts( indf, lfits, status );

      /* easiest to write a second piece of history information for header collation.
	 Stops having to worry about only getting a single HISTORY entry when NDF
	 wants to write a new entry every time the file is opened for UPDATE. */
      ndfHput("NORMAL",APPNAME, 1, 1, history,
	      0, 0, 0, indf, status );

      /* close file */
      ndfAnnul( &indf, status );

      /* free the copy and other objects */
      astAnnul( lfits );
      astAnnul( specwcs );

      /* we may want to set permissions on the file to stop unwary people overwriting it
	 or even the acquisition system itself. */
      if (*status == SAI__OK) {
	sysstat = chmod( fname, S_IRUSR | S_IWUSR| S_IRGRP | S_IROTH );
	if (sysstat == -1) {
	  *status = SAI__ERROR;
	  emsSyser( "ERRNO", errno );
	  emsSetc( "FILE", fname );
	  emsRep(" ","acsSpecCloseTS: Error setting permissions on file ^FILE: ^ERRNO",
		 status );
	}
      }

    }
  }

/* End the AST context. This annuls all AST Objects pointers created since
   the matching call to astBegin, except for any which have been exempted
   or exported. */
   astEnd;

  /* Reset AST status */
  astWatch( oldstat );

  return;
}

/********************************** Debug functions ********************/

#if SPW_DEBUG_LEVEL > 0
/* simply subtract two timeval structs and return the answer */
/* Does tp2 - tp1 */
static double duration ( struct timeval * tp1, struct timeval * tp2 ) {
  double diff = 0.0;
  diff = (tp2->tv_sec - tp1->tv_sec) +
    (tp2->tv_usec - tp1->tv_usec ) / 1E6;

  return diff;
}
#endif


AstFrameSet *specWcs( const AstFrameSet *fs, int ntime, const double times[], int * status ){

/*
*+
*  Name:
*     specWcs

*  Purpose:
*     Calculate frameset for spectrum time series.

*  Prototype:
*     AstFrameSet *specWcs( const AstFrameSet *fs, int ntime, const double times[],
*                int * status );

*  Description:
*     Returns a FrameSet in which the base Frame is a 3D GRID Frame, and
*     the current Frame has 3 axes in the order (spectrum,space,time).
*     The SpecFrame representing the spectral axis and its relationship
*     to GRID coords is read from the supplied FITS FrameSet. The time
*     axis is described using a MJD(TAI) TimeFrame, and its relationship
*     to GRID coords is specified by the supplied look-up table of time
*     values. The spatial axis is described by a simple 1D Frame with
*     Domain "SPACEINDEX" and is connected to the GRID coords via a
*     UnitMap.

*  Parameters:
*     fs = const AstFrameSet * (Given)
*        A pointer to the FrameSet read from the unmodified FITS headers
*        which define the spectral axis.
*     ntime = int (Given)
*        The number of time values supplied in "times".
*     times = const double [] (Given)
*        An array of "ntime" MJD values (in the TAI timescale), one for
*        each pixel along the time axis.
*     status = int * (Given & Returned)
*        Inherited status.

*  Returned Value:
*     specWcs = AstFrameSet *
*        3-D frameset.

*  Notes:
*     - The SpecFrame does not use the times in the time dimension

*  Authors:
*     DSB: David Berry (UCLan)
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     10-MAR-2006 (DSB):
*        Initial version (untested)
*     01-JUN-2006 (TIMJ):
*        Integrated into specwriter.

*-
*/

/* Local Variables: */
   AstCmpFrame *totfrm;
   AstCmpMap *totmap;
   AstFrame *axis, *specfrm, *spacefrm, *gridfrm;
   AstFrameSet *result;
   AstLutMap *timemap;
   AstMapping *specmap;
   AstTimeFrame *timefrm;
   AstUnitMap *spacemap;
   double tcopy[2];  /* local copy of time lut for when only 1 number present */
   const double *ltimes;  /* pointer to a time array */
   int nax, iax, iax_spec, ax_out[ NDF__MXDIM ];

/* Initialise. */
   result = NULL;

/* Check the global error status. */
   if( *status != SAI__OK ) return result;

/* Start an AST context so we do not need to annul AST pointers explicitly. */
   astBegin;

/* Check each axis of the current Frame in the FrameSet read from
   the FITS headers, looking for a spectral axis. We do this using
   astIsASpecFrame since this will pick up both SpecFrames and DSBSpecFrames
   (since a DSBSpecFrame "is a" SpecFrame). */
   specfrm = NULL;
   nax = astGetI( fs, "Naxes" );
   for( iax = 1; iax <= nax; iax++ ) {
      axis = astPickAxes( fs, 1, &iax, NULL );
      if( astIsASpecFrame( axis ) ) {
         specfrm = axis;
         iax_spec = iax;
         break;
      }
   }

/* Report an error if no spectral axis was found in the FITS header. */
   if( !specfrm ) {
      if( *status == SAI__OK ) {
         *status = SAI__ERROR;
         emsRep( "", "No spectral axis found in FITS header", status );
         goto L999;
      }
   }

/* We now assume that the spectral axis is connected to one and only one
   of the grid axes in the FITS header. We use astMapSplit to determine
   which grid axis this is, and to get the Mapping from the grid axis
   to the SpecFrame axis. We first invert the FrameSet (i.e. swap base
   and current Frames) since astMapSplit picks specified *inputs", but
   the SpecFrame is an "output" of the Mapping represented by the FrameSet.
   After inversion of the FrameSet, the SpecFrame will be one of the
   inputs, and can therefore be picked by astMapSplit. */
   astInvert( fs );
   astMapSplit( fs, 1, &iax_spec, ax_out, &specmap );

/* Invert the FrameSet again to return it to its original state */
   astInvert( fs );

/* Report an error if the assumption made above turned out not to be
   right. */
   if( !specmap ){
      if( *status == SAI__OK ) {
         *status = SAI__ERROR;
         emsRep( "", "The spectral axis depends on more than one pixel axis",
                 status );
         goto L999;
      }
   }

/* Invert "specmap" so that the forward transformation goes from grid
   coord to spectral coord. */
   astInvert( specmap );

/* We will use a simple Frame to describe the spatial axis, giving it the
   Domain name SPACEINDEX in order to distinguish it from the GRID Frame.
   The values on the spatial axis are just copies fo the grid coordinate,
   so create a UnitMap to connect the GRID Frame to the SPACEINDEX Frame. */
   spacefrm = astFrame( 1, "Domain=SPACEINDEX,Unit(1)=pixel,Label(1)=Receptor Number" );
   spacemap = astUnitMap( 1, "" );

/* We now have the SpecFrame, and the Mapping from grid coord to spectral
   coord. Now create a TimeFrame to describe MJD in the TAI timescale, and
   a LutMap which transforms grid coord into MJD (in days). The default
   TimeFrame attribute values give us what we want. Although we 
   can set the formatting to 2 decimal places (seconds) since we know that
   ACSIS will never read out faster than 50 ms.
*/
   timefrm = astTimeFrame( "" );
   if (ntime == 1) {
     /* a LutMap needs to numbers in its mapping so double up the
	first time if we only have one value. */
     tcopy[0] = times[0];
     tcopy[1] = times[1];
     ltimes = tcopy;
     ntime = 2;
   } else {
     ltimes = times;
   }
   timemap = astLutMap( ntime, ltimes, 1.0, 1.0, "" );

/* We now have the Frames and Mappings describing all the individual
   axes. Join all the Frames together into a CmpFrame (in the order spectral,
   spatial, time), and join all the Mappings together into a parallel
   CmpMap. */
   totfrm = astCmpFrame( astCmpFrame( specfrm, spacefrm, "" ), timefrm, "" );
   totmap = astCmpMap( astCmpMap( specmap, spacemap, 0, "" ), timemap, 0, "" );

/* Create a 3D GRID Frame. */
   gridfrm = astFrame( 3, "Domain=GRID,Title=FITS pixel coordinates" );
   astSet( gridfrm, "Unit(1)=pixel,Label(1)=FITS pixel axis 1" );
   astSet( gridfrm, "Unit(2)=pixel,Label(2)=FITS pixel axis 2" );
   astSet( gridfrm, "Unit(3)=pixel,Label(2)=FITS pixel axis 3" );

/* Create the FrameSet to return, initially containing just the above
   GRID Frame. */
   result = astFrameSet( gridfrm, "" );

/* Add the total Frame into the FrameSet using the total Mapping to
   connect it to the base (i.e. GRID) Frame. */
   astAddFrame( result, AST__BASE, totmap, totfrm );

/* Arrive here if an error occurs. Note, this is still inside the AST context
   delimited by astBegin/astEnd. */
L999:

/* If no error has occurred, export the resulting FrameSet pointer
   from the current AST context so that it will not be annulled by the
   following call to astEnd. If an error has occurred, annul it explicitly,
   in order to ensure we are returning a NULL pointer. */
   if( *status == SAI__OK ) {
      astExport( result );
   } else {
      result = astAnnul( result );
   }

/* End the AST context. This annuls all AST Objects pointers created since
   the matching call to astBegin, except for any which have been exempted
   or exported. */
   astEnd;

/* Return the resulting FrameSet. */
   return result;

}

/* inherited status is set to bad if the named file already exists */

static void checkNoFileExists( const char * file, int * status ) {

  struct stat fstat;     /* some where to store result from stat() */
  int err;               /* status code from stat() */

  if (*status != SAI__OK) return;

  /* There could be a race condition between the stat and subsequent
     but that should be impossible in normal DA usage so we ignore the
     possibility.
  */
  err = stat( file, &fstat);
  /* good error status means the file already exists */
  if (err == 0) {
    *status = SAI__ERROR;
    if (S_ISDIR(fstat.st_mode)) {
      emsSetc("T", "Directory");
    } else {
      emsSetc("T", "File");
    }
    emsSetc("F", file );
    emsRep(" ","^T ^F already exists. Not allowed to overwrite. "
	   "The data acquisition system is misconfigured - "
	   "please check your observation number is correct.",
	   status);
  }

}

/******************** kpgPtfts ***********************/

/* Stolen from kpgPtfts to minimize dependencies. Very bad. */



#define SZFITSCARD 80      /* Size of a FITS header card */
#define FITSSTR "80"       /* string representation of size of FITS */


/*
*+
*  Name:
*     kpgPtfts

*  Purpose:
*     Store FITS header information into an NDF

*  Language:
*     Starlink ANSI C

*  Invocation:
*     CALL KPG_PTFTS( INDF, FCHAN, STATUS )
*     kpgPtfts( int indf, AstFitsChan * fchan, int * status );

*  Description:
*     The routine stores the contents of an AST FitsChan into an
*     NDF by creating (or replacing) the FITS extension in the NDF.

*  Arguments:
*     indf = int (Given)
*        Identifier of NDF to receive the .FITS extension.
*     fchan = const AstFitsChan * (Given)
*        An AST pointer to a FitsChan which contains information about
*        the FITS header to be associated with the NDF.
*     status = int * (Given and Returned)
*        The global status.

*  Return Value:
*     Returns the status.

*  Notes:
*     - If a .MORE.FITS extension already exists it will be completely
*     replaced by this routine.

*  Copyright:
*     Copyright (C) 2005-2006 Particle Physics and Astronomy Research Council.
*     All Rights Reserved.

*  Licence:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     25-NOV-2005 (TIMJ):
*        Original version.
*     25-APR-2006 (TIMJ):
*        Finish.
*     {enter_changes_here}

*  Bugs:
*     {note_any_bugs_here}

*-
*/

static int kpgPtfts( int indf, const AstFitsChan * fchan, int * status ) {

  char card[SZFITSCARD+1];            /* A single FITS header card */
  int  fitsdim[1];  /* dimensions of FITS extension */
  HDSLoc * fitsloc = NULL;  /* Locator to FITS extension */
  char * fpntr;             /* Pointer to mapped FITS header */
  unsigned int i;           /* Loop counter */
  AstFitsChan * lchan; /* Local copy of FitsChan */
  unsigned int ncards;      /* Number of header cards */
  size_t nchars;            /* Actual size of FITS extension */
  int * oldstat;    /* Current status watched by AST */
  int result;               /* Result from astFindFits */
  int there = 0;    /* Is FITS extension there? */
  void * vpntr;     /* dummy void pointer */

  if ( *status != SAI__OK ) return *status;

  /* First need to look for a FITS extension */
  ndfXstat( indf, "FITS", &there, status );

  /* Remove it if it exists */
  if (there) {
    ndfXdel( indf, "FITS", status );
  }

  /* Make sure that we are checking AST status */
  oldstat = astWatch( status );

  /* Get local cloned copy of the FitsChan since we promised not
     to modify the supplied FitsChan */
  lchan = astCopy( fchan );

  /* Find out how many cards are present in the FitsChan */
  ncards = astGetI( lchan, "Ncard" );

  /* Rewind the FitsChan */
  astClear( lchan, "Card" );
    
  /* Create FITS extension */
  fitsdim[0] = ncards;
  ndfXnew(indf, "FITS", "_CHAR*" FITSSTR, 1, fitsdim, &fitsloc, status );

  /* Loop over all cards, inserting into extension -
     vpntr shenanigans fix strict-aliasing warning if casting
     &fpntr directly to void** */
  datMapV( fitsloc, "_CHAR*" FITSSTR, "WRITE", &vpntr, &nchars, status );
  fpntr = vpntr;

  if (*status == SAI__OK) {
    if ( ncards != nchars ) {
      *status = SAI__ERROR;
      emsSetu( "DM", nchars );
      emsSetu( "SZ",  ncards );
      emsRep("KPG_PTFTS_ERR",
	     "Bizarre error whereby number of cards in mapped FITS header (^DM) differs from number requested (^SZ)", status );
    }
  }

  if (*status == SAI__OK) {
    for (i = 1; i <= ncards; i++) {
      result = astFindFits( lchan, "%f", card, 1 );
      if (result) {
	strncpy( fpntr, card, SZFITSCARD );
	fpntr += SZFITSCARD;
      } else {
	break;
      }
    }
  }

  /* Cleanup */
  datUnmap( fitsloc, status );
  datAnnul( &fitsloc, status );

  /* Annul the local copy */
  astAnnul( lchan );

  /* Reset AST status */
  astWatch( oldstat );

  if (*status != SAI__OK)
    emsRep(" ", "kpgPtfts: Error writing FITS information to NDF", status );

  return *status;
}
