/*
*+
*  Name:
*     spwTest

*  Purpose:
*     Test the spectrum writing code.

*  Invocation:
*     make check

*  Language:
*     Starlink ANSI C

*  Description:
*     Simple test of spectrum writing.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     03-MAR-2006 (TIMJ):
*        Original version.
*     04-JAN-2006 (TIMJ):
*        Use jcmt/state interface
*     09-MAY-2007 (TIMJ):
*        OCSCONFIG added to open

*  Copyright:
*     Copyright (C) 2006,2007 Particle Physics and Astronomy Research Council.
*     Copyright (C) 2007 Science and Technology Facilities Council.
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

/* System includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>

#ifndef RAND_MAX
#define RAND_MAX INT_MAX
#endif

/* Starlink includes */
#include "ast.h"
#include "ndf.h"
#include "ems.h"
#include "prm_par.h"
#include "sae_par.h"
#include "star/hds.h"
#include "jcmt/state.h"

#define DEBUG_LEVEL 0

/* Local includes */
#include "specwrite.h"

/* Internal constants */
#define SPD 86400.0   /* seconds per day */

#define NSUBSYS 1     /* number of subsystems */
#define NCHAN   8192   /* number of channels */
#define NRECEP  16    /* Number of receptors */
#define DUMPRATE 20   /* in Hertz */


#define OBSLENGTH 1 /* Duration of the observation (seconds) */

/* Number of sequence steps in observation and the total number
   of spectra to write. SEQLEN is the length of a sequence in steps */
#define NSEQ  (int)( (float)OBSLENGTH * (float)DUMPRATE )
#define SEQLEN 10
#define NSPEC ( NSEQ * NRECEP * NSUBSYS )

/* Fraction of random lut to extract */
#define FRAC 0.5

/* data rate in MBps */
#define DATARATE ( DUMPRATE * NRECEP * NSUBSYS * NCHAN * 4 / (1024 * 1024) )

/* names of receptors */
static const char *recepnames[] = { "H01",
				    "H02",
				    "H03",
				    "H04",
				    "H05",
				    "H06",
				    "H07",
				    "H08",
				    "H09",
				    "H10",
				    "H11",
				    "H12",
				    "H13",
				    "H14",
				    "H15",
				    "H16",
};

static const char focal_station[] = "DIRECT";

static const char ocsconfig[] = "<OCS_CONFIG>\n</OCS_CONFIG>\n";

const float fplanex[] = { 1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.6f,8.2f,
			  9.2f,10.f,11.f,12.f,13.f,14.f,15.f,16.f };
const float fplaney[] = { -1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.6f,8.2f,
			  -9.2f,10.f,11.f,12.f,13.f,14.f,15.f,16.f };

/* Prototypes */
static int kpgGtfts( int indf, AstFitsChan ** fchan, int * status );
static double duration ( struct timeval * tp1, struct timeval * tp2, int * status );
static void writeSpectrum( const float * spectrum, unsigned int nsubsys,
			   JCMTState * record,
			   ACSISSpecHdr *spechdr,
			   size_t nchans[],
			   int feed, unsigned int * count, unsigned int lseq, double step_time,
			   int * status );

void MAIN_ ( void );
void MAIN_ ( void ) { };

int
main ( void ) {

  const unsigned int nsubsys = NSUBSYS;
  AstFitsChan * fitschan = NULL;
  AstFitsChan * fits[NSUBSYS];
  int status = SAI__OK;
  int i;
  float spectrum[NCHAN];
  size_t nchans[NSUBSYS];
  char lut[NRECEP*NSEQ];
  JCMTState record;
  ACSISSpecHdr spechdr;
  struct timeval tp1;
  struct timeval tp2;
  double diff;
  unsigned int c;
  unsigned int seq;
  int exstat;
  const double rts_step = 1.0 / DUMPRATE;
  int indf;
  int feed;

  emsBegin( &status );

  hdsTune("64BIT", 1, &status );

  /* Fill spectrum */
  for (i = 0; i < NCHAN; i++) {
    spectrum[i] = (float)i;
  }

  /* Create header - obtain from test file */
  ndfFind( NULL, "testhdr", &indf, &status );
  // ndfFind( NULL, "a20060920_00001_00_0001", &indf, &status );
  //ndfFind( NULL, "ac20061004_00071_01_01", &indf, &status );
  // ndfFind( NULL, "omc1_small", &indf, &status );
  kpgGtfts( indf, &fitschan, &status );
  ndfAnnul( &indf, &status );

  /* information */
  printf("Required data rate: %d MB/s\n", DATARATE );
  printf("Number of sequence steps: %d\n", NSEQ );
  printf("HDS TS cube. Writing %d spectra over %d seconds\n", (int)NSPEC, (int)OBSLENGTH);

  /* intialise the record */
  record.pol_ang = 3.14159;
  record.rts_num = 0;
  strcpy(record.rts_tasks, "SIMULATOR");
  strcpy(record.tcs_beam, " ");
  strcpy(record.smu_chop_phase, " ");
  spechdr.acs_tsys = 52.8;
  spechdr.acs_trx = 256.7;
  strcpy(record.tcs_source, "SCIENCE");
  *(record.tcs_tr_sys) = '\0';
  strcpy(record.acs_source_ro, "SPECTRUM_RESULT");
  record.acs_no_prev_ref = 1;
  record.acs_no_next_ref = 1;
  record.acs_no_ons = 1;
  record.acs_exposure = 1.0 / DUMPRATE;
  record.fe_lofreq = 245.6;
  record.fe_doppler = 1.0;
  spechdr.acs_feedx = 1.0;
  spechdr.acs_feedy = 2.0;
  record.jos_drcontrol = 2;
  

  /* Open NDF */
  for (i = 0; i < nsubsys; i++) {
    nchans[i] = NCHAN;
    fits[i] = fitschan;
  }

  acsSpecOpenTS( ".", 20060607, 53, NRECEP, nsubsys, recepnames,
		 focal_station, fplanex, fplaney, ocsconfig, &status);
  c = 0;
  spechdr.rts_endnum = 0;

  /* Initialise lut to zeroes */
  memset( lut, 0, NRECEP*NSEQ );

  /* go through and extract FRAC % of lut */
#if HAVE_SRANDOMDEV
  srandomdev();
#endif
  for (i=0; i < (int)(FRAC*NRECEP*NSEQ); i++) {
    if (status != SAI__OK) break;
    long rnd;
    while ( 1 ) {
      double frac = (double)random() / (double)RAND_MAX;
      double drnd = (double)(NSEQ * NRECEP) * frac;
      rnd = (int)drnd;
      if (lut[rnd] == 0 ) break; /* 0 means we have not sent it */
    }
    feed = rnd % NRECEP;
    seq = rnd / NRECEP;
#if DEBUG_LEVEL > 0
    printf("---------------Sending random seq %u feed %d\n", seq, feed);
#endif
    writeSpectrum( spectrum, nsubsys, &record, &spechdr, nchans, feed, &c, seq, rts_step, &status );
    lut[rnd]++;
  }

  /* now go through "lut" and send the remaining spectra */
  for (i=0; i< (int)(NSEQ*NRECEP); i++) {
    if (status != SAI__OK) break;
    if ( lut[i] == 0 ) {
      feed = i % NRECEP;
      seq = i / NRECEP;
#if DEBUG_LEVEL > 0
      printf("Sending seq %u feed %d\n", seq, feed);
#endif
      writeSpectrum( spectrum, nsubsys, &record, &spechdr, nchans, feed, &c, seq, rts_step, &status );
    }

  }

  /* Now close up */
  gettimeofday(&tp1, NULL);
  acsSpecCloseTS( fits, 0, &status );
  gettimeofday(&tp2, NULL);
  diff = duration( &tp1, &tp2, &status );
  printf("Time to close file = %.3f seconds\n", diff);

  hdsShow("LOCATORS", &status);
  hdsShow("FILES", &status);

  if (status == SAI__OK) {
    exstat = EXIT_SUCCESS;
  } else {
    printf("Status = %d\n", status);
    exstat = EXIT_FAILURE;
  }

  emsEnd( &status );

  return exstat;
}

static void writeSpectrum( const float spectrum[], unsigned int nsubsys, JCMTState *record,
			   ACSISSpecHdr * spechdr, size_t nchans[],
			   int feed, unsigned int * count, unsigned int seqnum, double step_time,
			   int * status ) {
  unsigned int j;
  float spectrum2[NCHAN];
  struct timeval tp1;
  struct timeval tp2;
  double diff;
  unsigned int nseqchunks;
  double ref_time = 53797.0; 
  double step_time_in_days;

  step_time_in_days = step_time / SPD;

  /* Set the feed and sequence number */
  spechdr->acs_feed = feed;
  record->rts_num = seqnum + 1;

  /* The end sequence number can be guessed */
  nseqchunks= seqnum / SEQLEN;
  spechdr->rts_endnum = (nseqchunks + 1) * SEQLEN;
  record->rts_end = step_time_in_days * seqnum + ref_time;
  record->tcs_tai = record->rts_end;

  /* tweak content */
  for (j=0; j < NCHAN; j++) {
    spectrum2[j] = spectrum[j] + (float)feed;
  }

  for (j = 0; j < nsubsys; j++) {
    gettimeofday(&tp1, NULL);
    (*count)++;
    /* printf("Writing spectrum sequence %u end %u\n",record->rts_num, record->rts_endnum); */
    acsSpecWriteTS(j+1, nchans[j], spectrum2, record, spechdr, status);
    gettimeofday(&tp2, NULL);
    diff = duration( &tp1, &tp2, status);
    if ( diff > 0.5 ) {
      printf("Scan %d  in feed %d of seq %u was written in %.3f seconds\n", *count, feed, seqnum, diff);
    }
    
  }

}

static double duration ( struct timeval * tp1, struct timeval * tp2, int * status ) {
  double diff = 0.0;
  if (*status != SAI__OK) return diff;

  diff = (tp2->tv_sec - tp1->tv_sec) +
    (tp2->tv_usec - tp1->tv_usec ) / 1E6;

  return diff;
}

/* Again - steal some code for reading the FITS header since it is in kaplibs
   and not in NDF */


#include <string.h>

#include "star/hds.h"
#include "ast.h"
#include "ndf.h"
#include "mers.h"
#include "sae_par.h"
#include "dat_par.h"
#include "ndf_err.h"
#include "dat_err.h"

#define SZFITSCARD 80      /* Size of a FITS header card */
#define FITSSTR "80"       /* string representation of size of FITS */

/*
*+
*  Name:
*     kpgGtfts

*  Purpose:
*     Obtain FITS header information from an NDF

*  Language:
*     Starlink ANSI C

*  Invocation:
*     CALL KPG_GTFTS( INDF, FCHAN, STATUS )
*     kpgGtfts( int indf, AstFitsChan ** fchan, int * status );

*  Description:
*     The routine reads the FITS extension from an NDF and returns an
*     AST pointer to a FitsChan which contains this information. The 
*     information may then be accessed using routines from the AST 
*     library (SUN/211).

*  Arguments:
*     indf = int (Given)
*        NDF identifier.
*     fchan = AstFitsChan ** (Returned)
*        An AST pointer to a FitsChan which contains information about
*        the FITS headers associated with the NDF.
*     status = int * (Given and Returned)
*        The global status.

*  Return Value:
*     Returns the status.

*  Notes:
*     - It is the caller's responsibility to annul the AST pointer
*     issued by this routine (e.g. by calling AST_ANNUL) when it is no
*     longer required.
*     - If this routine is called with STATUS set, then a value of
*     AST__NULL will be returned for the FCHAN argument, although no
*     further processing will occur. The same value will also be
*     returned if the routine should fail for any reason.
*     - Status is set to KPG__NOFTS if no FITS extension is found.

*  Copyright:
*     Copyright (C) 2005 Particle Physics and Astronomy Research Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*     
*     This program is distributed in the hope that it will be
*     useful,but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*     
*     You should have received a copy of the GNU General Public License
*     along with this program; if not, write to the Free Software
*     Foundation, Inc., 59 Temple Place,Suite 330, Boston, MA
*     02111-1307, USA

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     25-NOV-2005 (TIMJ):
*        Original version.
*     29-NOV-2005 (TIMJ):
*        Rename from ndfGtfts
*     {enter_changes_here}

*  Bugs:
*     {note_any_bugs_here}

*-
*/

static int kpgGtfts( int indf, AstFitsChan ** fchan, int * status ) {

  char   *card;               /* Pointer to start of current card */
  HDSLoc *fitsloc = NULL;     /* FITS HDS Locator in extension */
  hdsdim fitsdim[DAT__MXDIM]; /* Dimensionality of FITS extension */
  void   *fpntr = NULL;       /* Pointer to the mapped FITS header */
  unsigned int    i;          /* Loop counter */
  size_t ncards;              /* Number of header cards in extension */
  size_t nchars;              /* Number of characters in extension */
  int    ndim;                /* Number of dimensions in FITS array */
  int    *oldstat;            /* Current status watched by AST */
  int    there = 0;           /* Is FITS extension there? */
  char   type[DAT__SZTYP+1];  /* Data type of the FITS extension */

  /* make sure the fits chan is set to NULL on exit with bad status */
  *fchan = AST__NULL;

  if ( *status != SAI__OK ) return *status;

  /* First need to look for a FITS extension */
  ndfXstat( indf, "FITS", &there, status );

  if ( *status == SAI__OK ) {
    if (!there) {
      *status = SAI__ERROR;
      errRep( "KPG_GTFTS_NOF", "FITS extension is not present in NDF",
	      status );
    }
  }

  /* Get the locator to the FITS extension */
  ndfXloc( indf, "FITS", "READ", &fitsloc, status );

  /* Get the data type */
  datType( fitsloc, type, status );

  if ( *status == SAI__OK ) {
    if (strcmp(type, "_CHAR*" FITSSTR) != 0 ) {
      *status = DAT__TYPIN;
      msgSetc( "TYP", type );
      errRep( "KPG_GTFTS_TYP", "Data type of FITS extension is '^TYP' not '_CHAR*" FITSSTR "'", status );
    }
  }

  /* Determine the dimensionality of the FITS extension */
  datShape( fitsloc, DAT__MXDIM, fitsdim, &ndim, status );

  if ( *status == SAI__OK ) {
    if ( ndim != 1 ) {
      *status = DAT__DIMIN;
      msgSeti( "NDIM", ndim );
      errRep( "KPG_GTFTS_DIM", "Number of dimensions in FITS extension = ^NDIM but should be 1", status );
    }
  }

  /* Get number of FITS entries - should match fitsdim[0] */
  datSize( fitsloc, &ncards, status );

  if ( *status == SAI__OK ) {
    if ( ncards != (size_t)fitsdim[0] ) {
      *status = DAT__DIMIN;
      msgSeti( "DM", (int)fitsdim[0] );
      msgSeti( "SZ", (int)ncards );
      errRep( "KPG_GTFTS_SIZ","Bizarre error whereby the first dimension of the FITS extension (^DM) does not equal the size of the extension (^SZ)", status);
    }
  }

  /* Use datMapV to map the entire FITS array, then step through
     it 80 characters at a time until we have done all the cards.
     Note that there is no nul-terminator so we can not use
     astPutCards directly */
    
  datMapV( fitsloc, "_CHAR*" FITSSTR, "READ", &fpntr, &nchars, status );

  if ( *status == SAI__OK ) {
    if ( ncards != nchars ) {
      *status = DAT__DIMIN;
      msgSeti( "DM", (int)nchars);
      msgSeti( "SZ", (int)ncards );
      errRep( "KPG_GTFTS_SIZ2","Bizarre error whereby the number of elements mapped in the FITS extension (^DM) does not equal the size of the extension (^SZ)", status);
    }
  }

  /* Do not bother with AST stuff if status is bad */
  if ( *status == SAI__OK ) {

    /* Associate AST status with our status */
    oldstat = astWatch( status );

    /* Create a new FitsChan */
    *fchan = astFitsChan( NULL, NULL, "" );

    /* store pointer to start of string in new variable for iteration */
    card = fpntr;
  
    /* Extract headers 80 characters at a time. No nul-termination
       but astPutFits guarantees to only read 80 characters */
    for (i = 0; i < ncards; i++ ) {
      astPutFits( *fchan, card, 0 );
      card += SZFITSCARD;
    }

    /* Rewind the FitsChan */
    astClear( *fchan, "Card" );

    /* if status is bad annul the Fits Chan */
    if ( *status != SAI__OK ) astAnnul( *fchan );

    /* Reset AST status */
    astWatch( oldstat );

  }

  /* Clean up */
  datUnmap( fitsloc, status );
  datAnnul( &fitsloc, status );

  /* Report wrapper error message */
  if ( *status != SAI__OK ) {
    errRep( "KPG_GTFTS_ERR",
	    "KPG_GTFTS: Error obtaining FITS header from an NDF.",
	    status );
  }

  return *status;
}

