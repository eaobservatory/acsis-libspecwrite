
#define PACKAGE_VERSION  "V0.1-1"

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

/* Local includes */
#include "specwrite.h"

#define SPW_DEBUG 1

/* Largest file name allowed (including path) */
#define MAXFILE 1024


/* Global state variables */

/* Maximum number of subsystems we can handle 
   We know that ACSIS can have at most 4 spectral subsystems
   and they will not change during a single observation. */

#define MAXSUBSYS 4
const unsigned int maxsubsys = MAXSUBSYS;

/* Lookup table to enable us to name the files for each subsystem */
const char subsysprefix[MAXSUBSYS] = { 'a', 'b', 'c', 'd' };

/* Current sequence number (so that we know when we need to increment the time axis "counters")
   Note that the sequence number will not be sequential and will not start at zero. It will
   just guarantee to be at least the same as the previous value (for multi-receptor observations).
*/
static unsigned int curseq[MAXSUBSYS] = { 0,0,0,0 };

/* t- Position to write next spectrum - Fortran indexing */
static unsigned int counters[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Current size of NDF for each subsystem */
static unsigned int cursize[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Number of channels expected for each subsystem */
static unsigned int nchans_per_subsys[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Number of receptors in this particular observation */
static unsigned int nreceps_per_obs = 0;

/* Pointer into mapped data array for each subsystem */
static float * spectra[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* Array of Ast framesets dealing with the spectral axis. */

static AstFrameSet* frameset_cache[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* NDF identifier of spectrum file (if used). One NDF per subsystem. */
static int indf[MAXSUBSYS] = { NDF__NOID, NDF__NOID, NDF__NOID, NDF__NOID };

/* Information required to construct a filename */
static char datadir[MAXFILE+1];
unsigned int this_obs_num = 0;
unsigned int this_ut = 0; 
static unsigned int this_subscan = 0;


/* internal prototypes */
static void writeFitsChan( int indf, const AstFitsChan * fitschan, int * status );
static char * getFileName( const char * dir, unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status );
static char * getDirName( const char * dir, unsigned int yyyymmdd, 
			   unsigned int obsnum, int * status );
static char * getFileRoot( unsigned int yyyymmdd, unsigned int subsys,
			   unsigned int obsnum, unsigned int subscan, int * status );

static char * createSubScanDir( const char * dir, unsigned int yyyymmdd, 
				unsigned int obsnum,
				int * status );

static void
openNDF( unsigned int subsys, unsigned int nrecep, unsigned int nchans,
	 unsigned int nseq, const char * recepnames[], int * status );

static void
closeNDF( unsigned int subsys, int * status );

static void
resizeNDF( unsigned int subsys, unsigned int newsize, int * status );

static void
createExtensions( unsigned int subsys, unsigned int size, int * status );

static void
createACSISExtensions( unsigned int subsys, unsigned int size, unsigned int nrecep, 
		       const char *recepnames[], int * status );

static void resizeExtensions( unsigned int subsys, unsigned int newsize, int remap, 
			      int * status );
static void resizeACSISExtensions( unsigned int subsys, unsigned int newsize, int remap, 
			      int * status );
static void closeExtensions( unsigned int subsys, int * status );
static void closeACSISExtensions( unsigned int subsys, int * status );

static void writeRecord( unsigned int subsys, const ACSISRtsState * record, 
			 int * status );
static void writeRecepPos( unsigned int subsys, const ACSISRtsState * record, int * status );

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

/* Define the number of extensions we support */
#define NEXTENSIONS 39

/* Number of dimensions in output NDF */
#define NDIMS 3

/* definitions of dimensions */
#define CHANDIM 0
#define RECDIM  1
#define TDIM    2

/* Define indices for array of mapped pointers to extensions */
#define POL_ANG      0
#define RTS_NUM      1
#define RTS_STEP     2
#define RTS_END      3
#define RTS_TASKS    4
#define SMU_AZ_JIG_X 5
#define SMU_AZ_JIG_Y 6
#define SMU_X        7
#define SMU_Y        8
#define SMU_Z        9 
#define SMU_TR_JIG_X 10
#define SMU_TR_JIG_Y 11
#define TCS_AIRMASS  12
#define TCS_AZ_SYS   13
#define TCS_AZ_ANG   14
#define TCS_AZ_AC1   15
#define TCS_AZ_AC2   16
#define TCS_AZ_DC1   17
#define TCS_AZ_DC2   18
#define TCS_AZ_BC1   19
#define TCS_AZ_BC2   20
#define TCS_INDEX    21
#define TCS_SOURCE   22
#define TCS_TR_SYS   23
#define TCS_TR_ANG   24
#define TCS_TR_AC1   25
#define TCS_TR_AC2   26
#define TCS_TR_DC1   27
#define TCS_TR_DC2   28
#define TCS_TR_BC1   29
#define TCS_TR_BC2   30
#define ACS_SOURCE_RO    31
#define ACS_SOURCE_RP    32
#define ACS_DRCONTROL    33
#define ACS_TSYS         34
#define ACS_TRX          35
#define ENVIRO_REL_HUM   36
#define ENVIRO_PRESSURE  37
#define ENVIRO_AIR_TEMP  38

/* Definitions of HDS types associated with ACSISRtsStates struct. All these
   will be created in the file. */
static const char * hdsRecordNames[NEXTENSIONS][2] = 
  {
   { "_DOUBLE", "POL_ANG" },
   { "_INTEGER", "RTS_NUM" },
   { "_DOUBLE", "RTS_STEP" },
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
   { CHARTYP(SIZEOF_TCS_AZ_SYS), "TCS_AZ_SYS" },
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
   { CHARTYP(SIZEOF_ACS_SOURCE_RO), "ACS_SOURCE_RO" },
   { CHARTYP(SIZEOF_ACS_SOURCE_RP), "ACS_SOURCE_RP" },
   { "_INTEGER", "ACS_DRCONTROL" },
   { "_REAL",    "ACS_TSYS" },
   { "_REAL",     "ACS_TRX" },
   { "_REAL", "ENVIRO_REL_HUM" },
   { "_REAL", "ENVIRO_PRESSURE" },
   { "_REAL", "ENVIRO_AIR_TEMP" },
  };

/* Extension support */

/* Name of STATE and ACSIS extensions - some could almost be shared with SCUBA2... */
#define STATEEXT   "JCMTSTATE"
#define STATEEXTTYPE "RTS_ARR"

#define ACSISEXT   "ACSIS"
#define ACSISEXTTYP "ACSIS_COMP"

/* Locator to extensions */
static HDSLoc * extloc[MAXSUBSYS] = { NULL, NULL, NULL, NULL };
static HDSLoc * acsisloc[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* Flag to indicate that extensions are mapped */
static int extmapped[MAXSUBSYS] = { 0, 0, 0, 0 };
static int acsismapped[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Array of HDS locators to each of the extensions - different for each
 subsystem. */
static HDSLoc* extlocators[MAXSUBSYS][NEXTENSIONS];

/* Array of pointers to the mapped extensions */
static void * extdata[MAXSUBSYS][NEXTENSIONS];

/* HDS locator to RECEPPOS extension */
static HDSLoc* receppos_loc[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* pointer to mapped RECEPPOS extension */
static double * receppos_data[MAXSUBSYS] = { NULL, NULL, NULL, NULL };


/*********************** NDF "cube" FILE *************************************/

/* Number of sequences to increment file size by if it runs out of space */

/* May want this parameter to be settable as the number of spectra that
   we expect in a given time period */
#define MAXRECEP   16
#define MAXRATE    20
#define PRESIZETIME 10
#define NGROW  (MAXRATE * PRESIZETIME)

/* Number of bytes we should write before opening a new file */
#define MAXBYTES ( 100 * 1024 * 1024 )


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

*-
*/

void
acsSpecOpenTS( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
	       unsigned int nrecep, unsigned int nsubsys, 
	       const size_t nchans[], unsigned int nseq,
	       const char *recepnames[], int * status ) {

  unsigned int i;              /* Loop counter */
  char * sdir;                 /* subscan directory */

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
  for (i = 0; i < nsubsys; i++) {
    if (indf[i] != NDF__NOID) {
      *status = SAI__ERROR;
      emsSetu("I", i);
      emsRep("HDS_SPEC_OPENTS_ERR1",
	     "acsSpecOpenTS called, yet an NDF file is already open (subsystem ^I)",
	     status);
      return;
    }
  }

  /* Create the directory to receive each subscan */
  sdir = createSubScanDir( dir, yyyymmdd, obsnum, status );

  /* Store the global observation state */
  strncpy( datadir, sdir, MAXFILE );
  datadir[MAXFILE] = '\0';
  this_obs_num = obsnum;
  this_ut = yyyymmdd;

  /* Reset subscan number */
  this_subscan = 1;

  /* we need to store the receptor names somewhere locally */


  /* Need an NDF per subsystem */
  for (i = 0; i < nsubsys; i++) {

    openNDF( i, nrecep, nchans[i], nseq, recepnames, status );

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
acsSpecWriteTS( unsigned int subsys, const float spectrum[], 
		const ACSISRtsState* record,
		const AstFitsChan * freq, int * status ) {

  float * data; /* local copy of mapped pointer to spectrum */
  unsigned int offset;         /* offset into data array */
  int seqinc = 0;              /* did we increment sequence number? */
  unsigned int ngrow;          /* number of time slices to grow file */
  unsigned int reqnum;         /* number of time slices indicates by RTS sequence */
  AstSpecFrame * template;     /* Template of a specFrame */
  AstFrameSet  * foundframe;   /* Frameset from template to found frame */


  if (*status != SAI__OK) return;

  /* make sure that the subsys number is in range */
  if ( subsys >= maxsubsys ) {
    *status = SAI__ERROR;
    emsSetu("IN", subsys);
    emsSeti("MAX", maxsubsys-1);
    emsRep(" ","acsSpecWriteTS: Supplied subsystem number (^IN) exceeds max allowed (^MAX)", status);
    return;
  }

  /* Make sure the file is open */
  /* Check to see if we've already been called */
  if (indf[subsys] == NDF__NOID) {
    *status = SAI__ERROR;
    emsSetu("I", subsys);
    emsRep(" ",
	   "acsSpecWriteTS called, yet an NDF file has not been opened (subsystem ^I)",
	   status);
    return;
  }



  /* have we been given WCS? */
  if ( freq ) {
    /* Need to extract astrometry from the FitsChan */
  }

  /* if the sequence number has incremented we need to increase the t-axis counter */

  if ( record->rts_num > curseq[subsys] ) {
    /* store the new value */
    curseq[subsys] = record->rts_num;

    /* increment the counters value */
    counters[subsys]++;

    /* indicate that we did increment sequence number (and so can write header information) */
    seqinc = 1;

    /* See if we need to grow */
    /* We can grow either because we have suddenly realised we don't fit *or* because
       we have been told how many sequence steps to expect - calculate the required
       number to extend. (but we know at least 1) */
    reqnum = 1;

    if ( record->rts_endnum > record->rts_num ) {
      reqnum = record->rts_endnum - record->rts_num + 1;
    }

    if ( (counters[subsys] + reqnum - 1) > cursize[subsys] ) {
      /* resize NDF and all data arrays */

      /* work out how much to grow:
	 - use requested size if given and more than 1
	 - else use NGROW
	 - make sure we grow by at least counters[subsys]-cursize[subsys]
      */

      if (reqnum > 1) {
	ngrow = reqnum + counters[subsys] - cursize[subsys] - 1;
      } else {
	ngrow = NGROW;
      }

      /* Resize the NDF */
      resizeNDF( subsys, ngrow, status );
    }
  }

  /* copy in the data */
  if (*status == SAI__OK) {
    
    /* Calculate offset into array - number of spectra into the array times number of
       channels per spectrum. */
    offset = nchans_per_subsys[subsys] *
      (nreceps_per_obs * (counters[subsys] - 1) + record->acs_feed );

    data = spectra[subsys];
    memcpy( &(data[offset]), spectrum, nchans_per_subsys[subsys]*sizeof(float) );

    /* Store record data and receptor positions. Base record only updates each
       sequence step but recepot position should be written for all records. */
    if (seqinc) writeRecord( subsys, record, status );
    writeRecepPos( subsys, record, status );
    

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

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Loop over each NDF to write fits header and to close it */
  found = 0;
  for (i = 0; i < maxsubsys; i++) {
    if ( indf[i] != NDF__NOID ) {
      found = 1;
      closeNDF( i, status );
    }
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
    nchans_per_subsys[i] = 0;
  }

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
      flen = snprintf(rootname, MAXFILE, "a%c%d_%05d_%04d", subsysprefix[subsys], yyyymmdd, obsnum, subscan );
    } else {
      flen = snprintf(rootname, MAXFILE, "a%d_%05d", yyyymmdd, obsnum );
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
   - returns a pointer to static memory.
 */

static char * getDirName( const char * dir, unsigned int yyyymmdd,
			   unsigned int obsnum, int * status ) {

  static char dirname[MAXFILE]; /* buffer for dirname - will be returned */
  int flen;                        /* Length of string */
  char * root;                   /* Root filename */

  if (*status != SAI__OK) return NULL;
  root = getFileRoot( yyyymmdd, 0, obsnum, 0, status );
  if (*status != SAI__OK) return NULL;

  /* Form the file name - assume posix filesystem */
  flen = snprintf(dirname, MAXFILE, "%s/%s", dir, root );

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

/* Open an NDF file for that subsystem */

static void
openNDF( unsigned int subsys, unsigned int nrecep, unsigned int nchans,
	 unsigned int nseq, const char * recepnames[], int * status ) {

  void *datapntrs[] = { NULL };/* Array of mapped pointers for ndfMap */
  int itemp;                   /* Temp integer */
  int lbnd[NDIMS];             /* Lower pixel bounds */
  char *ndfname;               /* NDF filename */
  unsigned int ngrow;   /* Initial size to grow array */
  int place;                   /* NDF placeholder */
  int ubnd[NDIMS];             /* upper pixel bounds */
  char * history[1] = { "ACSIS Data Acquistion" };


#if SPW_DEBUG
  printf("Opening subsystem %d\n",subsys);
#endif
  /* Name the NDF component */
  ndfname = getFileName(datadir, this_ut, subsys, this_obs_num, this_subscan, status );

  /* Calculate bounds */
  ngrow = (nseq > 0 ? nseq : NGROW );
  lbnd[RECDIM] = 1;
  lbnd[CHANDIM] = 1;
  lbnd[TDIM] = 1;
  ubnd[RECDIM] = nrecep;
  ubnd[CHANDIM] = nchans;
  ubnd[TDIM] = ngrow;

#if SPW_DEBUG
  printf("Opening NDF file '%s' to default size of %u sequence steps\n", ndfname, ngrow);
#endif

  /* create the NDF */
  ndfPlace( NULL, ndfname, &place, status );
  ndfNew( "_REAL", NDIMS, lbnd, ubnd, &place, &(indf[subsys]), status );

  /* Update the cursize[] array and the nchans array */
  cursize[subsys] = ubnd[TDIM] - lbnd[TDIM] + 1;
  nchans_per_subsys[subsys] = ubnd[CHANDIM] - lbnd[CHANDIM] + 1;
  nreceps_per_obs = ubnd[RECDIM] - lbnd[RECDIM] + 1;

  /* History component */
  ndfHcre( indf[subsys], status );
  ndfHput("NORMAL","ACSIS-DA (" PACKAGE_VERSION ")", 1, 1, history,
	  0, 0, 0, indf[subsys], status );
  ndfHsmod( "DISABLED", indf[subsys], status );

  /* Map the data array */
  ndfMap(indf[subsys], "DATA", "_REAL", "WRITE", datapntrs, &itemp, status );

  /* Store the pointer */
  spectra[subsys] = datapntrs[0];

  /* Create the ACSIS extension that contains the receptor names and
     positions */
  createACSISExtensions( subsys, ngrow, nrecep, recepnames, status );

  /* Also need to create the header arrays and map those ! */
  createExtensions( subsys, ngrow, status );

  return;
}

/* 
   Create the extensions of specified size.
   Pointers stored in extdata.
   HDS Locators stored in extlocators.
*/

static void
createExtensions( unsigned int subsys, unsigned int size, int * status ) {

  int j;
  hdsdim dim[1];
  size_t ndim = 1;

  if (*status != SAI__OK) return;

  if (extmapped[subsys]) {
    *status = SAI__ERROR;
    emsRep(" ", "createExtensions: Attempting to create extensions that are already mapped", status );
    return;
  }

  /* Initial size */
  dim[0] = size;

  /* Create the extension */
  ndfXnew( indf[subsys], STATEEXT, STATEEXTTYPE, 0, NULL, &(extloc[subsys]), status ); 

  /* Loop and create. Can initialise HDS locator array safely */
  for (j=0; j < NEXTENSIONS; j++ ) {
    extlocators[subsys][j] = NULL;
    datNew( extloc[subsys], hdsRecordNames[j][1], hdsRecordNames[j][0],
	    ndim, dim, status );

    datFind( extloc[subsys], hdsRecordNames[j][1], &(extlocators[subsys][j]), status );

    datMap( extlocators[subsys][j], hdsRecordNames[j][0], "WRITE",
	    ndim, dim, &(extdata[subsys][j]), status );
    if ( *status != SAI__OK ) break;

  }

  if (*status == SAI__OK) extmapped[subsys] = 1;

  if (*status != SAI__OK)
    emsRep(" ", "Error creating JCMT state extension", status );

}

/*
  Resize the extensions to the supplied value.
  If remap is false, the arrays will not be remapped (so call at end to resize
  before annulling locators 
*/

static void
resizeExtensions( unsigned int subsys, unsigned int newsize, 
		  int remap, int * status ) {

  int j;
  hdsdim dim[1];
  size_t ndim = 1;

  if (*status != SAI__OK) return;

  dim[0] = newsize;

  /* Do all the unmapping. Then all the resizing then all the mapping */

  for (j=0; j < NEXTENSIONS; j++ ) {

    datUnmap( extlocators[subsys][j], status );
    if ( *status != SAI__OK ) break;
  }

  for (j=0; j < NEXTENSIONS; j++ ) {

    /* resize */
    datAlter( extlocators[subsys][j], 1, dim, status);
    if ( *status != SAI__OK ) break;
  }

  if (remap) {
    for (j=0; j < NEXTENSIONS; j++ ) {
      
      /* remap - assume this should be done after resizing all */
      datMap( extlocators[subsys][j], hdsRecordNames[j][0], "WRITE",
	      ndim, dim, &(extdata[subsys][j]), status );
      if ( *status != SAI__OK ) break;

    }
  }

  if (*status != SAI__OK)
    emsRep(" ", "Error resizing JCMT state extension", status );

}

/* Close down the extensions and free resources */

static void closeExtensions( unsigned int subsys, int * status ) {

  int j;

  if ( *status != SAI__OK ) return;

  if (counters[subsys] > 0) {
    resizeExtensions( subsys, counters[subsys], 0, status );
  }

  /* Free locators */
  for (j=0; j < NEXTENSIONS; j++) {
    datAnnul( &(extlocators[subsys][j]), status );
  }

  /* Close extension */
  datAnnul( &(extloc[subsys]), status );

  /* indicate that we are closed down */
  extmapped[subsys] = 0;

  /* delete the extension if we never wrote to it */
  if (counters[subsys] == 0) {
    ndfXdel(indf[subsys], STATEEXT,status);
  }

  if (*status != SAI__OK)
    emsRep(" ", "Error closing JCMT state extension", status );
}

/* Write ACSISRtsState to file */

static void writeRecord( unsigned int subsys, const ACSISRtsState * record,
			 int * status ) {

  unsigned int frame; /* position in data array */
  unsigned int offset;

  /* Can not think of anything clever to do */
  if ( *status != SAI__OK ) return;

  frame = counters[subsys] - 1;

  /* now copy */
  ((double *)extdata[subsys][POL_ANG])[frame] = record->pol_ang;
  ((int *)extdata[subsys][RTS_NUM])[frame] = record->rts_num;
  ((double *)extdata[subsys][RTS_STEP])[frame] = record->rts_step;
  ((double *)extdata[subsys][RTS_END])[frame] = record->rts_end;

  cnfExprt( record->rts_tasks,
	   (char *)extdata[subsys][RTS_TASKS]+
	    SIZEOF_RTS_TASKS*frame,
	    SIZEOF_RTS_TASKS );

  ((double *)extdata[subsys][SMU_AZ_JIG_X])[frame] = record->smu_az_jig_x;
  ((double *)extdata[subsys][SMU_AZ_JIG_Y])[frame] = record->smu_az_jig_y;
  ((double *)extdata[subsys][SMU_X])[frame] = record->smu_x;
  ((double *)extdata[subsys][SMU_Y])[frame] = record->smu_y;
  ((double *)extdata[subsys][SMU_Z])[frame] = record->smu_z;
  ((double *)extdata[subsys][SMU_TR_JIG_X])[frame] = record->smu_tr_jig_x;
  ((double *)extdata[subsys][SMU_TR_JIG_Y])[frame] = record->smu_tr_jig_y;
  ((double *)extdata[subsys][TCS_AIRMASS])[frame] = record->tcs_airmass;

  cnfExprt(record->tcs_az_sys,
	   (char *)extdata[subsys][TCS_AZ_SYS]+SIZEOF_TCS_AZ_SYS*frame, 
	   SIZEOF_TCS_AZ_SYS );

  ((double *)extdata[subsys][TCS_AZ_ANG])[frame] = record->tcs_az_ang;
  ((double *)extdata[subsys][TCS_AZ_AC1])[frame] = record->tcs_az_ac1;
  ((double *)extdata[subsys][TCS_AZ_AC2])[frame] = record->tcs_az_ac2;
  ((double *)extdata[subsys][TCS_AZ_DC1])[frame] = record->tcs_az_dc1;
  ((double *)extdata[subsys][TCS_AZ_DC2])[frame] = record->tcs_az_dc2;
  ((double *)extdata[subsys][TCS_AZ_BC1])[frame] = record->tcs_az_bc1;
  ((double *)extdata[subsys][TCS_AZ_BC2])[frame] = record->tcs_az_bc2;
  ((int *)extdata[subsys][TCS_INDEX])[frame] = record->tcs_index;
  cnfExprt ( record->tcs_source,
	     (char *)extdata[subsys][TCS_SOURCE]+SIZEOF_TCS_SOURCE*frame, 
	     SIZEOF_TCS_SOURCE );

  cnfExprt ( record->tcs_tr_sys,
	     (char *)extdata[subsys][TCS_TR_SYS]+SIZEOF_TCS_TR_SYS*frame, 
	     SIZEOF_TCS_TR_SYS );

  ((double *)extdata[subsys][TCS_TR_ANG])[frame] = record->tcs_tr_ang;
  ((double *)extdata[subsys][TCS_TR_AC1])[frame] = record->tcs_tr_ac1;
  ((double *)extdata[subsys][TCS_TR_AC2])[frame] = record->tcs_tr_ac2;
  ((double *)extdata[subsys][TCS_TR_DC1])[frame] = record->tcs_tr_dc1;
  ((double *)extdata[subsys][TCS_TR_DC2])[frame] = record->tcs_tr_dc2;
  ((double *)extdata[subsys][TCS_TR_BC1])[frame] = record->tcs_tr_bc1;
  ((double *)extdata[subsys][TCS_TR_BC2])[frame] = record->tcs_tr_bc2;

  cnfExprt( record->acs_source_ro,
	    (char*)extdata[subsys][ACS_SOURCE_RO]+ SIZEOF_ACS_SOURCE_RO*frame,
	    SIZEOF_ACS_SOURCE_RO );
  cnfExprt( record->acs_source_rp,
	    (char*)extdata[subsys][ACS_SOURCE_RP]+ SIZEOF_ACS_SOURCE_RP*frame,
	    SIZEOF_ACS_SOURCE_RP );

  ((int *)extdata[subsys][ACS_DRCONTROL])[frame] = record->acs_drcontrol;
  ((float *)extdata[subsys][ACS_TSYS])[frame] = record->acs_tsys;
  ((float *)extdata[subsys][ACS_TRX])[frame] = record->acs_trx;

  ((float *)extdata[subsys][ENVIRO_AIR_TEMP])[frame] = record->enviro_air_temp;
  ((float *)extdata[subsys][ENVIRO_PRESSURE])[frame] = record->enviro_pressure;
  ((float *)extdata[subsys][ENVIRO_REL_HUM])[frame] = record->enviro_rel_hum;

}

/* Create the .MORE.ACSIS extensions (that are not JCMT state structure members) */


static void
createACSISExtensions( unsigned int subsys, unsigned int size, unsigned int nrecep, 
		       const char *recepnames[], int * status ) {
  size_t maxlen = 0;
  size_t len;
  unsigned int i;
  char type[DAT__SZTYP+1];   /* constructed type string */
  HDSLoc * temploc = NULL;
  hdsdim dim[3];

  if (*status != SAI__OK) return;

  if (acsismapped[subsys]) {
    *status = SAI__ERROR;
    emsRep( " ", "createACSISExtensions: ACSIS extension already mapped. Can not create", status);
    return;
  }

  /* create ACSIS extension */
  ndfXnew( indf[subsys], ACSISEXT, ACSISEXTTYP, 0, NULL, &(acsisloc[subsys]), status );

  /* Need to create the following components:
     - RECEPTORS  _CHAR* array for each of the nrecep elements and their names. This fixed
       once written.
     - RECEPPOS   _DOUBLE (2 * size)   x and y positions (in tracking coordinates) for each
       receptor. This array grows in the same way as JCMTSTATE.
  */

  /* find the longest receptor name */
  if (recepnames != NULL) {
    for (i=0; i < nrecep; i++) {
      len = strlen( recepnames[i] );
      if (maxlen < len ) maxlen = len;
    }
  }

  /* Create the receptor component and store the names */
  if (maxlen > 0) {
    datCctyp( maxlen, type );
    dim[0] = nrecep;
    datNew( acsisloc[subsys], "RECEPTORS", type, 1, dim, status );
    datFind( acsisloc[subsys], "RECEPTORS", &temploc, status );
    if (recepnames != NULL) datPut1C( temploc, nrecep, recepnames, status );
    datAnnul( &temploc, status );
  }

  /* Now create the positions array and map it */
  dim[0] = 2;
  dim[1] = nrecep;
  dim[2] = size;
  datNew( acsisloc[subsys], "RECEPPOS", "_DOUBLE", 3, dim, status );
  datFind( acsisloc[subsys], "RECEPPOS", &(receppos_loc[subsys]), status );
  datMapD( receppos_loc[subsys], "WRITE", 3, dim, &(receppos_data[subsys]), status );

  if (*status != SAI__OK) {
    acsismapped[subsys] = 0;
    receppos_data[subsys] = NULL;
  } else {
    acsismapped[subsys] = 1;
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
resizeACSISExtensions( unsigned int subsys, unsigned int newsize, 
		       int remap, int * status ) {

  hdsdim dim[3];
  size_t ndim = 3;
  int actdim;

  if (*status != SAI__OK) return;

  dim[0] = newsize;

  if (acsismapped[subsys])
    datUnmap( receppos_loc[subsys], status );

  /* Get the current bounds */
  datShape( receppos_loc[subsys], ndim, dim, &actdim, status );

  if (*status != SAI__OK && ndim != (size_t)actdim) {
    *status = SAI__ERROR;
    emsSeti( "AD", actdim);
    emsSetu( "ND", (unsigned long) ndim );
    emsRep(" ", "Dims mismatch in ACSIS extension. ^AD != ^ND", status);
  }

  /* resize */
  dim[ndim-1] = newsize;
  datAlter( receppos_loc[subsys], ndim, dim, status);

  if (remap) {
    /* remap - assume this should be done after resizing all */
    datMapD( receppos_loc[subsys], "WRITE",
	     ndim, dim, &(receppos_data[subsys]), status );
  }

  if (*status != SAI__OK)
    emsRep(" ", "Error resizing ACSIS extension", status );

}

/* Close down the ACSIS extension and free resources */

static void closeACSISExtensions( unsigned int subsys, int * status ) {

  if ( *status != SAI__OK ) return;

  if (counters[subsys] > 0) {
    resizeACSISExtensions( subsys, counters[subsys], 0, status );
  }

  /* Free locators */
  datAnnul( &(receppos_loc[subsys]), status );

  /* delete the receptor positions if never written */
  if (counters[subsys] == 0) {
    datErase(acsisloc[subsys], "RECEPPOS", status );
  }

  /* Close extension */
  datAnnul( &(acsisloc[subsys]), status );

  acsismapped[subsys] = 0;

  if (*status != SAI__OK)
    emsRep(" ", "Error closing ACSIS extension", status );

}

/* Write coordinate positions to ACSIS extension */
static void writeRecepPos( unsigned int subsys, const ACSISRtsState * record, int * status ) {
  unsigned int offset;
  double *posdata;

  if (*status != SAI__OK) return;

  /* Calculate offset into data array */
  offset = 2 * ( (nreceps_per_obs * (counters[subsys]-1) ) + record->acs_feed );

  posdata = receppos_data[subsys];
  if (posdata != NULL) {
    posdata[offset] = record->acs_feedx;
    posdata[offset+1] = record->acs_feedy;
  } else {
    *status = SAI__ERROR;
    emsRep( " ", "Attempted to write receptor positions but no data array mapped",
	    status );
  }

}

/* Close NDF for a particular subsystem */

static void
closeNDF( unsigned int subsys, int * status ) {

  int itemp;                /* Temp integer */
  int lbnd[NDIMS];          /* Lower bounds of NDF */
  int ubnd[NDIMS];          /* upper bounds of NDF */

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Unmap */
#if SPW_DEBUG
  printf("Unmap current NDF final time\n");
#endif
  TIMEME("Final unmap", ndfUnmap( indf[subsys], "DATA", status ););

  /* Shrink file to actual size */
  ndfBound(indf[subsys], NDIMS, lbnd, ubnd, &itemp, status );
  if (counters[subsys] > 0) {
    ubnd[TDIM] = lbnd[TDIM] + counters[subsys] - 1;
  } else {
    ubnd[TDIM] = lbnd[TDIM];
  }

#if SPW_DEBUG
  printf("Setting final bounds. Resize to %lld spectra\n", (unsigned long long)ubnd[TDIM]);
#endif
  TIMEME( "Final set bounds", ndfSbnd(NDIMS, lbnd, ubnd, indf[subsys], status ););

  /* Close extensions */
  TIMEME( "Final extension resize", closeExtensions( subsys, status ););
  closeACSISExtensions( subsys, status );

  /* Close file */
  ndfAnnul( &(indf[subsys]), status );

#if SPW_DEBUG
  printf("Wrote %d sequence steps to subsystem %d (max was %d)\n", counters[subsys], subsys,
	 cursize[subsys]);
#endif

  /* Force globals to be reset */
  cursize[subsys] = 0;
  counters[subsys] = 0;
  spectra[subsys] = NULL;
  indf[subsys] = NDF__NOID;

  if (*status != SAI__OK) {
    emsSetu("S", subsys);
    emsRep( " ", "Error closing Spectrum NDF file subsystem ^S", status );
  }

}


/* Resize a specific NDF */

static void
resizeNDF( unsigned int subsys, unsigned int newsize, int * status ) {

  int ubnd[NDIMS];
  int lbnd[NDIMS];
  unsigned int nchans;
  unsigned int nreceps;
  unsigned int newt;
  int itemp;
  void *datapntrs[] = { NULL };/* Array of mapped pointers for ndfMap */

  /* Unmap the data array */
#if SPW_DEBUG
  printf("Unmap data array in preparation for resize\n");
#endif
  TIMEME( "ndfUnmap", ndfUnmap(indf[subsys], "DATA", status ););

  /* Get the existing bounds */
  ndfBound(indf[subsys], NDIMS, lbnd, ubnd, &itemp, status );

  if (*status == SAI__OK && itemp != NDIMS) {
    *status = SAI__ERROR;
    emsSeti("N", itemp);
    emsSeti("ND", NDIMS);
    emsRep(" ", "acsSpecWriteTS: Bizarre internal error. Ndims is ^N not ^ND",
	   status);
  }
    
  nchans = ubnd[CHANDIM] - lbnd[CHANDIM] + 1;
  if (*status == SAI__OK && nchans != nchans_per_subsys[subsys]) {
    *status = SAI__ERROR;
    emsSetu("UB", nchans);
    emsSetu("NC", nchans_per_subsys[subsys] );
    emsRep(" ", "acsSpecWriteTS: Bizzare internal error. Nchans is ^UB not ^NC",
	   status);
  }

  nreceps = ubnd[RECDIM] - lbnd[RECDIM] + 1;
  if (*status == SAI__OK && nreceps != nreceps_per_obs) {
    *status = SAI__ERROR;
    emsSetu("UB", nreceps);
    emsSetu("NR", nreceps_per_obs);
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
  TIMEME("Set NDF bounds", ndfSbnd( NDIMS, lbnd, ubnd, indf[subsys], status ););

  /* map data array again */
#if SPW_DEBUG
  printf("Remap the data array\n");
#endif
  TIMEME("Remp NDF", ndfMap( indf[subsys], "DATA", "_REAL", "WRITE", datapntrs, &itemp, status ););
  spectra[subsys] = datapntrs[0];

  /* Resize the extensions */
  TIMEME("Resize extensions", resizeExtensions( subsys, newt, 1, status  ););
  TIMEME("Resize coords", resizeACSISExtensions( subsys, newt, 1, status  ););
  
  /* Update cursize */
  cursize[subsys] = newt;

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
