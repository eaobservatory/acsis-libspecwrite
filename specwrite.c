
#define PACKAGE_VERSION  "V0.1-1"

/* System includes */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/* Starlink includes */
#include "sae_par.h"
#include "star/hds.h"
#include "ndf.h"
#include "ast.h"
#include "ems.h"
#include "mers.h"
#include "prm_par.h"

/* Local includes */
#include "hdsspec.h"

/* Global state variables */

/* Locator of the root HDS container file */

static HDSLoc * locator = NULL;

/* Maximum number of subsystems we can handle 
   We know that ACSIS can have at most 4 spectral subsystems
   and they will not change during a single observation. */

#define MAXSUBSYS 4
const unsigned int maxsubsys = MAXSUBSYS;

/* Current spectrum counter - increment *after* writing a new spectrum.
   Scalar is used for HDS .I1, .I2 writing.
   Array is used for multi-subsystem time series NDFs.
*/
static unsigned int counter = 1;

/* Position to write next spectrum - C indexing */
static unsigned int counters[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Current size of NDF for each subsystem */
static int cursize[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Number of channels expected for each subsystem */
static unsigned int nchans_per_subsys[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Pointer into mapped data array for each subsystem */
static float * spectra[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* Array of Ast framesets dealing with the spectral axis. */

static AstFrameSet* frameset_cache[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* NDF identifier of spectrum file (if used). One NDF per subsystem. */
static int indf[MAXSUBSYS] = { NDF__NOID, NDF__NOID, NDF__NOID, NDF__NOID };


/* internal prototypes */
static void writeFitsChan( int indf, const AstFitsChan * fitschan, int * status );
static char * getFileName( const char * dir, unsigned int yyyymmdd, 
			 unsigned int obsnum, int * status );
static void openHDSContainer( const char * dir, unsigned int yyyymmdd,
			      unsigned int obsnum, int * status );
static void
createExtensions( unsigned int subsys, unsigned int size, int * status );

static void resizeExtensions( unsigned int subsys, unsigned int newsize, int remap, 
			      int * status );
static void closeExtensions( unsigned int subsys, int * status );
static void writeRecord( unsigned int subsys, const ACSISRecord * record, 
			 int * status );

/* Largest file name allowed (including path) */
#define MAXFILE 1024

/* Function to put quotes around a symbol so that we can do
   CPP string concatenation */
#define myxstr(s) mystr(s)
#define mystr(s) #s
#define CHARTYP(s) "_CHAR*" myxstr(s)

/* Define the number of extensions we support */
#define NEXTENSIONS 44

/* Define indices for array of mapped pointers to extensions */
#define POL_ANG      0
#define RTS_NUM      1
#define RTS_STEP     2
#define RTS_END      3
#define RTS_TASKS    4
#define SMU_AZ_OFF_X 5
#define SMU_AZ_OFF_Y 6
#define SMU_X        7
#define SMU_Y        8
#define SMU_Z        9 
#define SMU_TR_OFF_X 10
#define SMU_TR_OFF_Y 11
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
#define ACS_RECEPTOR     31
#define ACS_SOURCE_RO    32
#define ACS_SOURCE_RP    33
#define ACS_DRCONTROL    34
#define ACS_TSYS         35
#define ACS_TRX          36
#define WVM_TH       37
#define WVM_T12      38
#define WVM_T42      39
#define WVM_T78      40
#define WVM_TW       41
#define WVM_QUAL     42
#define WVM_TIME     43  

/* Definitions of HDS types associated with ACSISRecords struct. All these
   will be created in the file. */
static const char * hdsRecordNames[NEXTENSIONS][2] = 
  {
   { "_DOUBLE", "POL_ANG" },
   { "_INTEGER", "RTS_NUM" },
   { "_DOUBLE", "RTS_STEP" },
   { "_DOUBLE", "RTS_END" },
   { CHARTYP(SIZEOF_RTS_TASKS), "RTS_TASKS" },
   { "_DOUBLE", "SMU_AZ_OFF_X" },
   { "_DOUBLE", "SMU_AZ_OFF_Y" },
   { "_DOUBLE", "SMU_X" },
   { "_DOUBLE", "SMU_Y" },
   { "_DOUBLE", "SMU_Z" },
   { "_DOUBLE", "SMU_TR_OFF_X" },
   { "_DOUBLE", "SMU_TR_OFF_Y" },
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
   { CHARTYP(SIZEOF_ACS_RECEPTOR), "ACS_RECEPTOR" },
   { CHARTYP(SIZEOF_ACS_SOURCE_RO), "ACS_SOURCE_RO" },
   { "_INTEGER", "ACS_SOURCE_RP" },
   { "_INTEGER", "ACS_DRCONTROL" },
   { "_REAL",    "ACS_TSYS" },
   { "_REAL",     "ACS_TRX" },
   { "_REAL", "WVM_TH" },
   { "_REAL", "WVM_T12" },
   { "_REAL", "WVM_T42" },
   { "_REAL", "WVM_T78" },
   { "_REAL", "WVM_TW" },
   { "_INTEGER", "WVM_QUAL" },
   { "_REAL", "WVM_TIME" }
  };

/* Extension support */

/* Name of extension - could almost be shared with SCUBA2... */
static char extname[] = "JCMTSTATE";
static char exttype[] = "RTS_ARR";

/* Locator to extension */
static HDSLoc * extloc[MAXSUBSYS] = { NULL, NULL, NULL, NULL };

/* Flag to indicate that extensions are mapped */
static int extmapped[MAXSUBSYS] = { 0, 0, 0, 0 };

/* Array of HDS locators to each of the extensions - different for each
 subsystem. */
static HDSLoc* extlocators[MAXSUBSYS][NEXTENSIONS];

/* Array of pointers to the mapped extensions */
static void * extdata[MAXSUBSYS][NEXTENSIONS];

/*********************** NDF "cube" FILE *************************************/

/* Number of spectra to increment file size by if it runs out of space */

/* May want this parameter to be settable as the number of spectra that
   we expect in a given time period */
#define MAXRECEP   16
#define MAXRATE    20
#define PRESIZETIME 10
#define NGROW  (MAXRECEP * MAXRATE * PRESIZETIME)


/*
*+
*  Name:
*     hdsSpecOpenTS

*  Purpose:
*     Open NDF file for writing time series spectra

*  Invocation:
*     hdsSpecOpenTS( const char * dir, unsigned int yyyymmdd, 
*                    unsigned int obsnum, unsigned int nrecep,
*                    unsigned int nsubsys, const size_t nchans[], 
*                    int * status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to prepare the file for output. It must
*     be called before calling hdsSpecWriteTS. It will be pre-sized to receive
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
*     a file is open. Call hdsSpecClose to close the file.
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
hdsSpecOpenTS( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
	       unsigned int nrecep, unsigned int nsubsys, 
	       const size_t nchans[], 
	       const char *recepnames[], int * status ) {

  void *datapntrs[] = { NULL };/* Array of mapped pointers for ndfMap */
  unsigned int i;              /* Loop counter */
  int itemp;                   /* Temporary integer */
  int lbnd[2];                 /* Lower pixel bounds */
  char ndfname[DAT__SZNAM+1];  /* NDF filename */
  int nlen;                    /* Number of characters in NDF component name */
  int place;                   /* NDF placeholder */
  int ubnd[2];                 /* upper pixel bounds */

  char * history[1] = { "ACSIS Data Acquistion" };

  /* Return immediately if status is bad */
  if (*status != SAI__OK) return;

  /* Validate the subsystem count */
  if (nsubsys > maxsubsys) {
    *status = SAI__ERROR;
    emsSeti("NIN", nsubsys);
    emsSeti("MAX", maxsubsys);
    errRep("HDS_SPEC_OPENTS_ERR0",
	   "hdsSpecOpenTS: number of subsystems supplied (^NIN) exceeds expected maximum of ^MAX", status);
    return;
  }

  /* Check to see if we've already been called */
  for (i = 0; i < nsubsys; i++) {
    if (indf[i] != NDF__NOID) {
      *status = SAI__ERROR;
      emsSeti("I", i);
      errRep("HDS_SPEC_OPENTS_ERR1",
	     "hdsSpecOpenTS called, yet an NDF file is already open (subsystem ^I)",
	     status);
      return;
    }
  }

  /* Open the container file */
  openHDSContainer( dir, yyyymmdd, obsnum, status );

  /* Need an NDF per subsystem */
  for (i = 0; i < nsubsys; i++) {

    printf("Opening subsystem %d\n",i);

    /* Name the NDF component */
    nlen = snprintf(ndfname, DAT__SZNAM+1, "SUBSYS%u", i );

    if (nlen > DAT__SZNAM) {
      *status = SAI__ERROR;
      errRep("HDS_SPEC_OPENTS_ERR2",
	     "Buffer overflow when forming NDF component name", status);
      return;
    }

    /* Calculate bounds */
    lbnd[0] = 1;
    lbnd[1] = 1;
    ubnd[0] = nchans[i];
    ubnd[1] = NGROW;

    /* create the NDF */
    ndfPlace( locator, ndfname, &place, status );
    ndfNew( "_REAL", 2, lbnd, ubnd, &place, &(indf[i]), status );

    /* Update the cursize[] array and the nchans array */
    cursize[i] = ubnd[1];
    nchans_per_subsys[i] = nchans[i];

    /* History component */
    ndfHcre( indf[i], status );
    ndfHput("NORMAL","ACSIS-DA (" PACKAGE_VERSION ")", 1, 1, history,
	    0, 0, 0, indf[i], status );
    ndfHsmod( "DISABLED", indf[i], status );

    /* Map the data array */
    ndfMap(indf[i], "DATA", "_REAL", "WRITE", datapntrs, &itemp, status );

    /* Store the pointer */
    spectra[i] = datapntrs[0];

    /* Also need to create the header arrays and map those ! */
    createExtensions( i, NGROW, status );

  }

  /* Complete */
  return;

}

/*
*+
*  Name:
*     hdsSpecWriteTS

*  Purpose:
*     Write a spectrum to an HDS time-series file.

*  Invocation:
*     hdsSpecWriteTS( unsigned int subsys, const float spectrum[], 
*                     const ACSISRecord * record, const ACSISFreqInfo * freq,
*                     int *status);

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to write a spectrum to an HDS container
*     that had previously been opened using hdsSpecOpen.

*  Arguments:
*     subsys = unsigned int (Given)
*        Subsystem used for this spectrum (start counting at 0).
*        Much match the nchans[] array given to hdsSpecOpenTS.
*     spectrum = float[nchan] (Given)
*        Spectrum itself.
*     record = const ACSISRecord * (Given)
*        Header information associated with this spectrum.
*     freq = const ACSISFreqInfo * (Given)
*        Frequency information associated with this spectrum.
*        Only required for the first spectrum in any subsystem.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.

*  Notes:
*     - Must have previously called hdsSpecOpenTS.

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
hdsSpecWriteTS( unsigned int subsys, const float spectrum[], 
		const ACSISRecord* record,
		const ACSISFreqInfo * freq, int * status ) {

  float * data; /* local copy of mapped pointer to spectrum */
  void *datapntrs[] = { NULL };/* Array of mapped pointers for ndfMap */
  int itemp;                   /* Temporary integer */
  int lbnd[2];                 /* Lower pixel bounds */
  int ubnd[2];                 /* upper pixel bounds */
  unsigned int offset;         /* offset into data array */

  if (*status != SAI__OK) return;

  /* make sure that the subsys number is in range */
  if ( subsys >= maxsubsys ) {
    *status = SAI__ERROR;
    emsSeti("IN", subsys);
    emsSeti("MAX", maxsubsys-1);
    errRep(" ","hdsSpecWriteTS: Supplied subsystem number (^IN) exceeds max allowed (^MAX)", status);
    return;
  }

  /* Make sure the file is open */
  /* Check to see if we've already been called */
  if (indf[subsys] == NDF__NOID) {
    *status = SAI__ERROR;
    emsSeti("I", subsys);
    errRep(" ",
	   "hdsSpecWriteTS called, yet an NDF file has not been opened (subsystem ^I)",
	   status);
    return;
  }



  /* write .WCS? */
  if (counters[subsys] == 0) {
  }

  /* See if we need to grow */
  if (counters[subsys] > cursize[subsys]-1) {
    /* resize NDF and all data arrays */

    /* Unmap the data array */
    ndfUnmap(indf[subsys], "DATA", status );

    /* Get the existing bounds */
    ndfBound(indf[subsys], 2, lbnd, ubnd, &itemp, status );

    if (*status == SAI__OK && itemp != 2) {
      *status = SAI__ERROR;
      emsSeti("N", itemp);
      errRep(" ", "hdsSpecWriteTS: Bizarre internal error. Ndims is ^N not 2",
	     status);
    }
    
    if (*status == SAI__OK && ubnd[0] != nchans_per_subsys[subsys]) {
      *status = SAI__ERROR;
      emsSeti("UB", ubnd[0]);
      emsSeti("NC", nchans_per_subsys[subsys] );
      errRep(" ", "hdsSpecWriteTS: Bizzare internal error. Nchans is ^UB not ^NC",
	     status);
    }

    /* increment */
    ubnd[1] += NGROW;

    /* set new bounds */
    ndfSbnd( 2, lbnd, ubnd, indf[subsys], status );

    /* map data array again */
    ndfMap( indf[subsys], "DATA", "_REAL", "WRITE", datapntrs, &itemp, status );
    spectra[subsys] = datapntrs[0];

    /* Resize the extensions */
    resizeExtensions( subsys, ubnd[1], 1, status  );

    /* Update cursize */
    cursize[subsys] = ubnd[1];
  }

  /* copy in the data */
  if (*status == SAI__OK) {
    
    /* Calculate offset into array */
    offset = nchans_per_subsys[subsys] * counters[subsys];

    data = spectra[subsys];
    memcpy( &(data[offset]), spectrum, nchans_per_subsys[subsys]*sizeof(float) );

    /* Store record data */
    writeRecord( subsys, record, status );

    /* increment the position */
    counters[subsys]++;

  }

  return;
}

/*
*+
*  Name:
*     hdsSpecCloseTS

*  Purpose:
*     Write FITS header and close HDS file.

*  Invocation:
*     hdsSpecCloseTS( const AstFitsChan * fits[], int *status );

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
*     - Must have previously called hdsSpecOpenTS.
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
hdsSpecCloseTS( const AstFitsChan * fits[], int * status ) {

  unsigned int i;           /* Loop counter */
  int found = 0;            /* Found an open NDF? */
  int itemp;                /* Temp integer */
  int lbnd[2];              /* Lower bounds of NDF */
  int ubnd[2];              /* upper bounds of NDF */

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Is the file opened? */
  if (locator == NULL) {
    *status = SAI__ERROR;
    errRep( " ",
	    "hdsSpecCloseTS: No file is open. Has hdsSpecOpen been called?", status );
    return;
  }

  /* Loop over each NDF to write fits header and to close it */
  found = 0;
  for (i = 0; i < maxsubsys; i++) {
    if ( indf[i] != NDF__NOID ) {
      found = 1;

      /* Unmap */
      ndfUnmap( indf[i], "DATA", status );

      /* Shrink file to actual size */
      ndfBound(indf[i], 2, lbnd, ubnd, &itemp, status );
      ubnd[1] = counters[i];
      ndfSbnd(2, lbnd, ubnd, indf[i], status );

      /* Close extensions */
      closeExtensions( i, status );

      /* FITS header */
      writeFitsChan( indf[i], fits[i], status );

      /* Close file */
      ndfAnnul( &(indf[i]), status );

      printf("Wrote %d spectra to subsystem %d (max was %d)\n", counters[i], i,
	     cursize[i]);

    }
  }

  /* report error if not found any open NDFs */
  if (*status == SAI__OK && !found) {
    *status = SAI__ERROR;
    errRep(" ", "hdsSpecCloseTS: Failed to find open NDF components", status );
  }

  /* Close the file */
  hdsClose( &locator, status);



  /* Force globals to be reset */
  for (i = 0; i < maxsubsys; i++) {
    cursize[i] = 0;
    counters[i] = 0;
    spectra[i] = NULL;
    indf[i] = NDF__NOID;
    nchans_per_subsys[i] = 0;
  }


}

/*********************** HDS CONTAINER FILE *************************************/

/*
*+
*  Name:
*     hdsSpecOpen

*  Purpose:
*     Open HDS container file and prepare for writing of spectra

*  Invocation:
*     hdsSpecOpen( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
*                  int * status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to prepare the file for output. It must
*     be called before calling hdsSpecWrite.

*  Arguments:
*     dir = const char * (Given)
*        Directory to write the file.
*     yyyymmdd = unsigned int (Given)
*        UT date in YYYYMMDD format. Used to construct the file name.
*     obsnum = unsigned int (Given)
*        Current observation number. Used to construct the file name.
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
*     a file is open. Call hdsSpecClose to close the file.
*     - The file created by this routine will be of the form
*       aYYYYMMDD_NNNNN.sdf where NNNNN is the zero padded observation number.

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
hdsSpecOpen( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
	     int * status ) {

  /* Return immediately if status is bad */
  if (*status != SAI__OK) return;

  /* Open the HDS container file */
  openHDSContainer( dir, yyyymmdd, obsnum, status );

  /* Complete */
  return;

}

/*
*+
*  Name:
*     hdsSpecWrite

*  Purpose:
*     Write a spectrum to an HDS file.

*  Invocation:
*     hdsSpecWrite( size_t nchan, const float spectrum[], const ACSISRecord * record,
*                   const ACSISFreqInfo * freq, int *status);

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to write a spectrum to an HDS container
*     that had previously been opened using hdsSpecOpen.

*  Arguments:
*     nchan = size_t (Given)
*        Number of channels in the spectrum.
*     spectrum = float[nchan] (Given)
*        Spectrum itself.
*     record = const ACSISRecord * (Given)
*        Header information associated with this spectrum.
*     freq = const ACSISFreqInfo * (Given)
*        Frequency information associated with this spectrum.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.

*  Notes:
*     - Must have previously called hdsSpecOpen.

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
hdsSpecWrite( size_t nchan, const float spectrum[], const ACSISRecord* record,
	      const ACSISFreqInfo * freq, int * status ) {

  void * datapntrs[] = { NULL }; /* Array of pointers from ndfMap */
  float * data;     /* _REAL pointer to file data array */
  int indf= NDF__NOID; /* NDF identifier */
  int lbnd[3];      /* Lower bounds of "cube" */
  int nelem;        /* Number of elements in output file */
  char ndfname[20]; /* Name of NDF spectrum "Innn" */
  int place;        /* Placeholder for NDF */
  int ubnd[3];      /* Upper bounds of "cube" */

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Is the file opened? */
  if (locator == NULL) {
    *status = SAI__ERROR;
    errRep( "HDS_SPEC_WRITE_ERR0",
	    "No file is open. Has hdsSpecOpen been called?", status );
    return;
  }

  /* Work out the name we are going to use for the spectrum */
  sprintf(ndfname, "I%d", counter );

  /* This spectrum will be written as a 1x1xN cube to simplify mosaicking */
  lbnd[0] = 1;
  ubnd[0] = 1;
  lbnd[1] = 1;
  ubnd[1] = 1;
  lbnd[2] = 1;     /* Should we put central channel at 0? */
  ubnd[2] = nchan;

  /* Create the NDF */
  ndfPlace( locator, ndfname, &place, status );
  ndfNew( "_REAL", 3, lbnd, ubnd, &place, &indf, status );

  /* Get a pointer to the data */
  ndfMap( indf, "DATA", "_REAL", "WRITE", datapntrs, &nelem, status );
  data = datapntrs[0];

  /* Validate */
  if (*status == SAI__OK && nelem != nchan ) {
    *status = SAI__ERROR;
    emsSeti("NC", nchan );
    emsSeti("NE", nelem );
    emsRep( "HDS_WRITE_SPEC_ERR2",
	    "Number of elements in file (^NE) differs from number supplied (^NC)",
	    status );
  }

  /* Copy the data into the file */
  if (*status == SAI__OK) memcpy( data, spectrum, nchan * sizeof(float) );

  /* Close the file */
  ndfAnnul( &indf, status );

  /* Okay - increment counter */
  if (*status == SAI__OK) counter++;

  return;
}

/*
*+
*  Name:
*     hdsSpecClose

*  Purpose:
*     Write FITS header and close HDS file.

*  Invocation:
*     hdsSpecClose( const AstFitsChan * fits, int *status );

*  Language:
*     Starlink ANSI C

*  Description:
*     This function must be used to close the file after all spectra have been
*     written. The FITS header is written.

*  Arguments:
*     fits = const AstFitsChan * (Given)
*        FITS header.
*     status = int * (Given & Returned)
*        Inherited status.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
*        Original version.

*  Notes:
*     - Must have previously called hdsSpecOpen.

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
hdsSpecClose( const AstFitsChan * fits, int * status ) {

  unsigned int i;           /* Loop counter */
  int itemp;                /* Temp integer */
  int indf;                 /* HEADER NDF identifier */
  const int lbnd[] = { 1 }; /* Lower bound of HEADER NDF */
  int place;                /* Place holder for NDF in HDS system */
  void * datapntrs[] = { NULL }; /* Array of pointers from ndfMap */
  unsigned char * tempdata;          /* UBYTE pointer */
  const int ubnd[] = { 1 }; /* Upper bound of HEADER NDF */

  /* Always check status on entry */
  if (*status != SAI__OK) return;

  /* Is the file opened? */
  if (locator == NULL) {
    *status = SAI__ERROR;
    errRep( "HDS_SPEC_WRITE_ERR0",
	    "No file is open. Has hdsSpecOpen been called?", status );
    return;
  }

  /* Write the FITS header */
  if (fits) {
    /* Create new NDF 1x1 called "HEADER" */
    ndfPlace( locator, "HEADER", &place, status );
    ndfNew( "_UBYTE", 1, lbnd, ubnd, &place, &indf, status );

    /* Map it to force it to be defined */
    ndfMap(indf, "DATA", "_UBYTE", "WRITE", datapntrs, &itemp, status );
    tempdata = datapntrs[0];
    *tempdata = VAL__BADUB;

    /* Write the FitsChan to the NDF */
    writeFitsChan( indf, fits, status );

    /* Close the files */
    ndfAnnul( &indf, status );

  }

  /* Close the file */
  hdsClose( &locator, status);

  /* Force globals to be reset */
  counter = 1;

  /* Annul spec frame sets */
  for (i = 0; i < maxsubsys; i++ ) {
    if ( frameset_cache[i] != NULL ) 
      frameset_cache[i] = astAnnul( frameset_cache[i] );
  }

}

/*********************** HELPER FUNCTIONS (PRIVATE) *****************************/

/* Internal function to write a FitsChan to an NDF - taken from
   kpgPtfts */

static void writeFitsChan( int indf, const AstFitsChan * fits, int * status ) {

  char card[81];            /* A single FITS header card */
  HDSLoc * fitsloc = NULL;  /* Locator to FITS extension */
  char * fpntr;             /* Pointer to mapped FITS header */
  int i;                    /* Loop counter */
  int ncards;               /* Number of header cards */
  size_t nchars;            /* Actual size of FITS extension */
  int result;               /* Result from astFindFits */

  if (*status != SAI__OK) return;

  /* Find out how many header cards we have */
  ncards = astGetI( fits, "Ncard" );

  /* Rewind the FitsChan */
  astClear( fits, "Card" );
    
  /* Create FITS extension */
  ndfXnew(indf, "FITS", "_CHAR*80", 1, &ncards, &fitsloc, status );

  /* Loop over all cards, inserting into extension */
  datMapV( fitsloc, "_CHAR*80", "WRITE", (void**)&fpntr, &nchars, status );

  if (*status == SAI__OK) {
    if ( ncards != nchars ) {
      *status = SAI__ERROR;
      emsSeti( "DM", nchars );
      emsSeti( "SZ",  ncards );
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

}

/* Form the file name
   - returns a pointer to static memory.
 */

static char * getFileName( const char * dir, unsigned int yyyymmdd,
			   unsigned int obsnum, int * status ) {

  static char filename[MAXFILE]; /* buffer for filename - will be returned */
  int flen;                        /* Length of string */

  /* Form the file name - assume posix filesystem */
  flen = snprintf(filename, MAXFILE, "%s/a%d_%05d", dir, yyyymmdd, obsnum );

  if (flen >= MAXFILE) {
    *status = SAI__ERROR;
    emsSeti("SZ", MAXFILE );
    emsSeti("N", obsnum );
    emsSeti("UT", yyyymmdd );
    errRep("HDS_SPEC_OPEN_ERR1",
	   "Error forming filename. Exceeded buffer size of ^SZ chars for scan ^N on UT ^UT", status );
    return NULL;
  }

  return filename;
}

/* Open a root HDS container file of the correct name.
   The locator itself is stored in the global "locator" variable.
*/

static void
openHDSContainer( const char * dir, unsigned int yyyymmdd, unsigned int obsnum,
	     int * status ) {

  const hdsdim dims[] = { 1 }; /* HDS dimensions (always ignored) */ 
  char *filename;              /* HDS root filename */

  /* Return immediately if status is bad */
  if (*status != SAI__OK) return;

  /* Check to see if we've already been called */
  if (locator != NULL) {
    *status = SAI__ERROR;
    errRep("HDS_SPEC_OPEN_ERR0",
	   "hdsSpecOpen called, yet an HDS file is already open", status);
    return;
  }

  /* Form the file name - assume posix filesystem */
  filename = getFileName( dir, yyyymmdd, obsnum, status );

  /* Try to open the file */
  hdsNew( filename, "ACSIS_SPEC", "CONTAINER", 0, dims, &locator, status );

  if (*status != SAI__OK) {
    locator = NULL;
    return;
  }

  /* Complete */
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
    errRep(" ", "createExtensions: Attempting to create extensions that are already mapped", status );
    return;
  }

  /* Initial size */
  dim[0] = size;

  /* Create the extension */
  ndfXnew( indf[subsys], extname, exttype, 0, NULL, &(extloc[subsys]), status ); 

  /* Loop and create. Can initialise HDS locator array safely */
  for (j=0; j < NEXTENSIONS; j++ ) {
    extlocators[subsys][j] = NULL;
    datNew( extloc[subsys], hdsRecordNames[j][1], hdsRecordNames[j][0],
	    ndim, dim, status );
    if ( *status != SAI__OK ) break;

    datFind( extloc[subsys], hdsRecordNames[j][1], &(extlocators[subsys][j]), status );
    if ( *status != SAI__OK ) break;

    datMap( extlocators[subsys][j], hdsRecordNames[j][0], "WRITE",
	    ndim, dim, &(extdata[subsys][j]), status );
    if ( *status != SAI__OK ) break;

  }

  if (*status == SAI__OK) extmapped[subsys] = 1;

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

}

/* Close down the extensions and free resources */

static void closeExtensions( unsigned int subsys, int * status ) {

  int j;

  if ( *status != SAI__OK ) return;

  resizeExtensions( subsys, counters[subsys], 0, status );

  /* Free locators */
  for (j=0; j < NEXTENSIONS; j++) {
    datAnnul( &(extlocators[subsys][j]), status );
  }

  /* Close extension */
  datAnnul( &(extloc[subsys]), status );

  /* indicate that we are closed down */
  if (*status == SAI__OK) extmapped[subsys] = 0;

}

/* Write ACSISRecord to file */

static void writeRecord( unsigned int subsys, const ACSISRecord * record, 
			 int * status ) {

  unsigned int frame; /* position in data array */

  /* Can not think of anything clever to do */
  if ( *status != SAI__OK ) return;

  frame = counters[subsys];

  /* now copy */
  ((double *)extdata[subsys][POL_ANG])[frame] = record->pol_ang;
  ((int *)extdata[subsys][RTS_NUM])[frame] = record->rts_num;
  ((double *)extdata[subsys][RTS_STEP])[frame] = record->rts_step;
  ((double *)extdata[subsys][RTS_END])[frame] = record->rts_end;

  cnfExprt( record->rts_tasks,
	   (char *)extdata[subsys][RTS_TASKS]+
	    SIZEOF_RTS_TASKS*frame,
	    SIZEOF_RTS_TASKS );

  ((double *)extdata[subsys][SMU_AZ_OFF_X])[frame] = record->smu_az_off_x;
  ((double *)extdata[subsys][SMU_AZ_OFF_Y])[frame] = record->smu_az_off_y;
  ((double *)extdata[subsys][SMU_X])[frame] = record->smu_x;
  ((double *)extdata[subsys][SMU_Y])[frame] = record->smu_y;
  ((double *)extdata[subsys][SMU_Z])[frame] = record->smu_z;
  ((double *)extdata[subsys][SMU_TR_OFF_X])[frame] = record->smu_tr_off_x;
  ((double *)extdata[subsys][SMU_TR_OFF_Y])[frame] = record->smu_tr_off_y;
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

  cnfExprt( record->acs_receptor,
	    (char*)extdata[subsys][ACS_RECEPTOR]+ SIZEOF_ACS_RECEPTOR*frame,
	    SIZEOF_ACS_RECEPTOR );

  cnfExprt( record->acs_source_ro,
	    (char*)extdata[subsys][ACS_SOURCE_RO]+ SIZEOF_ACS_SOURCE_RO*frame,
	    SIZEOF_ACS_SOURCE_RO );

  ((int *)extdata[subsys][ACS_SOURCE_RP])[frame] = record->acs_source_rp;
  ((int *)extdata[subsys][ACS_DRCONTROL])[frame] = record->acs_drcontrol;
  ((float *)extdata[subsys][ACS_TSYS])[frame] = record->acs_tsys;
  ((float *)extdata[subsys][ACS_TRX])[frame] = record->acs_trx;

  ((float *)extdata[subsys][WVM_TH])[frame] = record->wvm_th;
  ((float *)extdata[subsys][WVM_T12])[frame] = record->wvm_t12;
  ((float *)extdata[subsys][WVM_T42])[frame] = record->wvm_t42;
  ((float *)extdata[subsys][WVM_T78])[frame] = record->wvm_t78;
  ((float *)extdata[subsys][WVM_TW])[frame] = record->wvm_tw;
  ((int *)extdata[subsys][WVM_QUAL])[frame] = record->wvm_qual;
  ((float *)extdata[subsys][WVM_TIME])[frame] = record->wvm_time;

}


/*******************************************************************************/

/* Test MAIN */

/* 
gcc -g -O2 -I/star/include -I. hdsspec.c -L/star/lib `ndf_link` `one_link` -L/usr/local/g95/lib/gcc-lib/powerpc-apple-darwin8.2.0/4.0.1/ -lf95
*/

#define SPD 86400.0   /* seconds per day */

#define NSUBSYS 2
#define NCHAN  4096
#define NRECEP  16
/* in Hertz */
#define DUMPRATE 20

/* Duration of the observation (seconds) */
#define OBSLENGTH 60

#define NSPEC ( OBSLENGTH * NRECEP * DUMPRATE )

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

  hdsSpecOpenTS( ".", 20060607, 53, NRECEP, nsubsys, nchans, NULL, &status );
  for (i = 0; i < NSPEC; i++) {
    /* increment the sequence number every NRECEP spectra */
    if (i%NRECEP == 0) {
      record.rts_num ++;
      record.rts_end += step_time_in_days;
    }
    gettimeofday(&tp1, NULL);
    for (j = 0; j < nsubsys; j++) {
      hdsSpecWriteTS(j, spectrum, &record, NULL, &status);
    }
    gettimeofday(&tp2, NULL);
    diff = (tp2.tv_sec - tp1.tv_sec) +
      (tp2.tv_usec - tp1.tv_usec ) / 1E6;
    if ( diff > 0.5 ) {
      printf("Scan %d was written in %.3f seconds\n", i, diff);
    }
  }
  hdsSpecCloseTS( fits, &status );

  hdsShow("LOCATORS", &status);
  hdsShow("FILES", &status);

  return EXIT_SUCCESS;

  /* Open file */
  printf("HDS.In\n");
  hdsSpecOpen( ".", 20060607, 52, &status );

  /* Write a spectrum */
  for (i=0; i < NSPEC; i++) {
    hdsSpecWrite( NCHAN, spectrum, NULL, NULL, &status );
  }

  /* Close and write header */
  hdsSpecClose( fitschan, &status );

  /* End */
  return EXIT_SUCCESS;
}

void MAIN_ ( void ) { };
