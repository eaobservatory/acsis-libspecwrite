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

/* System includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/* Starlink includes */
#include "ast.h"
#include "ems.h"
#include "prm_par.h"
#include "sae_par.h"
#include "star/hds.h"


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
#define NSEQ  ( OBSLENGTH * DUMPRATE )
#define SEQLEN 10
#define NSPEC ( NSEQ * NRECEP * NSUBSYS )

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

static double duration ( struct timeval * tp1, struct timeval * tp2, int * status );

void MAIN_ ( void ) { };

int
main ( void ) {

  const unsigned int nsubsys = NSUBSYS;
  AstFitsChan * fitschan = NULL;
  AstFitsChan * fits[NSUBSYS];
  int status = SAI__OK;
  int i,j;
  float spectrum[NCHAN];
  float spectrum2[NCHAN];
  size_t nchans[NSUBSYS];
  ACSISRtsState record;
  double step_time_in_days;
  struct timeval tp1;
  struct timeval tp2;
  double diff;
  unsigned int c;
  unsigned int seq;
  int exstat;
  unsigned int lseq;

  emsBegin( &status );

  hdsTune("64BIT", 1, &status );

  /* Fill spectrum */
  for (i = 0; i < NCHAN; i++) {
    spectrum[i] = (float)i;
  }

  /* Create header */
  fitschan = astFitsChan( NULL, NULL, "" );
  astSetFitsS( fitschan, "INSTRUME", "ACSIS", "Instrument", 0 );
  astSetFitsS( fitschan, "TELESCOP", "JCMT", "Telescope", 0 );
  astShow( fitschan );

  printf("Required data rate: %d MB/s\n", DATARATE );
  printf("Number of sequence steps: %d\n", NSEQ );
  printf("HDS TS cube. Writing %d spectra over %d seconds\n", (int)NSPEC, (int)OBSLENGTH);

  /* intialise the record */
  record.pol_ang = 3.14159;
  record.rts_num = 0;
  record.rts_step = 1.0 / DUMPRATE;
  step_time_in_days = record.rts_step / SPD;
  record.rts_end  = 53797.0;
  strcpy(record.rts_tasks, "SIMULATOR");
  record.acs_trx = VAL__BADR;
  record.acs_tsys = VAL__BADR;

  /* Open NDF */
  for (i = 0; i < nsubsys; i++) {
    nchans[i] = NCHAN;
    fits[i] = fitschan;
  }

  acsSpecOpenTS( ".", 20060607, 53, NRECEP, nsubsys, nchans, (NSEQ/10), recepnames, &status);
  c = 0;
  record.rts_endnum = 0;
  for (seq = 0; seq < NSEQ; seq++) {
    lseq = seq;

    /* Increment sequence number in record */
    record.rts_num = lseq + 1;

    /* assume new sequence every 10th step */
    if ( lseq % SEQLEN == 0) record.rts_endnum += SEQLEN;

    record.rts_end += step_time_in_days;
    for (i = 0; i < NRECEP; i++) {
      record.acs_feed = i;
      record.acs_feedx = 1.0;
      record.acs_feedy = VAL__BADD;

      /* tweak content */
      for (j=0; j < NCHAN; j++) {
	spectrum2[j] = spectrum[j] + (float)i;
      }

      for (j = 0; j < nsubsys; j++) {
	gettimeofday(&tp1, NULL);
	c++;
	acsSpecWriteTS(j, spectrum2, &record, NULL, &status);
	gettimeofday(&tp2, NULL);
	diff = duration( &tp1, &tp2, &status);
	if ( diff > 0.5 ) {
	  printf("Scan %d  in feed %d of seq %u was written in %.3f seconds\n", c, i, lseq, diff);
	}

      }

    }
  }
  gettimeofday(&tp1, NULL);
  acsSpecCloseTS( fits, &status );
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



static double duration ( struct timeval * tp1, struct timeval * tp2, int * status ) {
  double diff = 0.0;
  if (*status != SAI__OK) return diff;

  diff = (tp2->tv_sec - tp1->tv_sec) +
    (tp2->tv_usec - tp1->tv_usec ) / 1E6;

  return diff;
}
