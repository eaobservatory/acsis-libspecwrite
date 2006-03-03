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
#include "prm_par.h"
#include "sae_par.h"
#include "star/hds.h"


/* Local includes */
#include "specwrite.h"

/* Internal constants */
#define SPD 86400.0   /* seconds per day */

#define NSUBSYS 2
#define NCHAN  4096
#define NRECEP  16
/* in Hertz */
#define DUMPRATE 20

/* Duration of the observation (seconds) */
#define OBSLENGTH 60

#define NSPEC ( OBSLENGTH * NRECEP * DUMPRATE )

void MAIN_ ( void ) { };

int
main ( void ) {

  const unsigned int nsubsys = NSUBSYS;
  AstFitsChan * fitschan = NULL;
  AstFitsChan * fits[NSUBSYS];
  int status = SAI__OK;
  int i,j;
  float spectrum[NCHAN];
  size_t nchans[NSUBSYS];
  ACSISRecord record;
  double step_time_in_days;
  struct timeval tp1;
  struct timeval tp2;
  double diff;

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

  acsSpecOpenTS( ".", 20060607, 53, NRECEP, nsubsys, nchans, NULL, &status );
  for (i = 0; i < NSPEC; i++) {
    /* increment the sequence number every NRECEP spectra */
    if (i%NRECEP == 0) {
      record.rts_num ++;
      record.rts_end += step_time_in_days;
    }
    gettimeofday(&tp1, NULL);
    for (j = 0; j < nsubsys; j++) {
      acsSpecWriteTS(j, spectrum, &record, NULL, &status);
    }
    gettimeofday(&tp2, NULL);
    diff = (tp2.tv_sec - tp1.tv_sec) +
      (tp2.tv_usec - tp1.tv_usec ) / 1E6;
    if ( diff > 0.5 ) {
      printf("Scan %d was written in %.3f seconds\n", i, diff);
    }
  }
  acsSpecCloseTS( fits, &status );

  hdsShow("LOCATORS", &status);
  hdsShow("FILES", &status);

  return EXIT_SUCCESS;
}


