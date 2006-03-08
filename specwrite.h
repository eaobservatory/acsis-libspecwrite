
#ifndef INCLUDED_HDS_SPEC_H
#define INCLUDED_HDS_SPEC_H

/*
*+
*  Name:
*     hdsspec.h

*  Purpose:
*     Define the public C interface to the HDS spectrum writing ACSIS functions

*  Invocation:
*     #include "hdsspec.h"

*  Language:
*     C Include file

*  Description:
*     This module defines the public interface to the HDS spectrum
*     writing functions used by ACSIS.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     27-FEB-2006 (TIMJ):
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

/* Some sizes used for dimensioning the string arrays */
#define SIZEOF_RTS_TASKS  80
#define SIZEOF_RTS_ERRS   80
#define SIZEOF_TCS_AZ_SYS 16
#define SIZEOF_TCS_SOURCE 32
#define SIZEOF_TCS_TR_SYS 16
#define SIZEOF_ACS_RECEPTOR   10
#define SIZEOF_ACS_SOURCE_RO  32

/* Public data structures */

/* State structures - based on sc2head from SCUBA-2 software*/
typedef struct ACSISRtsState {
  double pol_ang;
  unsigned int rts_num;
  double rts_step;
  double rts_end;    /* MJD TAI of end of sequence step */
  char   rts_tasks[SIZEOF_RTS_TASKS + 1];
  char   rts_errs[SIZEOF_RTS_ERRS + 1];
  double smu_az_off_x;
  double smu_az_off_y;
  double smu_x;
  double smu_y;
  double smu_z;
  double smu_tr_off_x;
  double smu_tr_off_y;
  double tcs_airmass;
  char   tcs_az_sys[SIZEOF_TCS_AZ_SYS+1];
  double tcs_az_ang;
  double tcs_az_ac1;
  double tcs_az_ac2;
  double tcs_az_dc1;
  double tcs_az_dc2;
  double tcs_az_bc1;
  double tcs_az_bc2;
  int tcs_index;
  char tcs_source[SIZEOF_TCS_SOURCE+1];
  char tcs_tr_sys[SIZEOF_TCS_TR_SYS+1];
  double tcs_tr_ang;
  double tcs_tr_ac1;
  double tcs_tr_ac2;
  double tcs_tr_dc1;
  double tcs_tr_dc2;
  double tcs_tr_bc1;
  double tcs_tr_bc2;
  float acs_tsys;
  float acs_trx;
  float wvm_th;
  float wvm_t12;
  float wvm_t42;
  float wvm_t78;
  float wvm_tw;
  int wvm_qual;
  float wvm_time;
  int   acs_feed;
  char  acs_source_ro[SIZEOF_ACS_SOURCE_RO+1];
  int   acs_source_rp;
  int   acs_spec_window_id;
  int   acs_drcontrol;
} ACSISRtsState;

/* Spectral scale */
typedef struct ACSISFreqInfo {
  double iffreq;
  double finc;
  double fcen;
  double refchan;
} ACSISFreqInfo;

/* NDF versions of the Spectrum writing */

/* Initialise the file for writing */
void acsSpecOpenTS( const char * dir,
		    unsigned int yyyymmdd,
		    unsigned int obsnum,
		    unsigned int nrecep,
		    unsigned int nsubsys,
		    const size_t nchans[],
		    unsigned int nseq,
		    const char* recepnames[],
		    int * status );

/* Write a spectrum to the file */
void acsSpecWriteTS( unsigned int subsys,
		     const float spectrum[],
		     const ACSISRtsState * state,
		     const ACSISFreqInfo * freq,
		     int * status );

/* Close the file */
void acsSpecCloseTS( const AstFitsChan * fits[],
		   int * status );



/* INCLUDE_HDS_SPEC_H */
#endif
