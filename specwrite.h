
#ifndef INCLUDED_HDS_SPEC_H
#define INCLUDED_HDS_SPEC_H

/*
*+
*  Name:
*     hdsspec.h

*  Purpose:
*     Define the public C interface to the HDS spectrum writing ACSIS functions

*  Invocation:
*     #include "acsis/specwrite.h"

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
*     21-APR-2006 (TIMJ):
*        nchans now used in specWrite not specOpen.

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
#define SIZEOF_TCS_SOURCE 32
#define SIZEOF_TCS_TR_SYS 16
#define SIZEOF_ACS_RECEPTOR   5
#define SIZEOF_ACS_SOURCE_RO  16

/* Public data structures */

/* State structures - based on sc2head from SCUBA-2 software*/
typedef struct ACSISRtsState {
  double pol_ang;
  unsigned int rts_num;
  unsigned int rts_endnum;  /* Highest number expected in this sequence */
  double rts_end;    /* MJD TAI of end of sequence step */
  char   rts_tasks[SIZEOF_RTS_TASKS + 1];
  char   rts_errs[SIZEOF_RTS_ERRS + 1];
  int    smu_jig_index;
  double smu_az_jig_x;
  double smu_az_jig_y;
  double smu_az_chop_x;
  double smu_az_chop_y;
  double smu_x;
  double smu_y;
  double smu_z;
  double smu_tr_jig_x;
  double smu_tr_jig_y;
  double smu_tr_chop_x;
  double smu_tr_chop_y;
  double tcs_airmass;
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
  float enviro_rel_hum;
  float enviro_pressure;
  float enviro_air_temp;
  unsigned int acs_feed;  /* Feed number */
  float acs_tsys;
  float acs_trx;
  char  acs_source_ro[SIZEOF_ACS_SOURCE_RO+1];
  int   acs_spec_window_id;
  int   acs_drcontrol;
  int   acs_no_prev_refs;
  int   acs_no_next_refs;
  int   acs_no_ons;
  float acs_exposure;
  double acs_feedx;  /* Y coordinate of feed "acs_feed" */
  double acs_feedy;  /* Y coordinate of feed "acs_feed" */
} ACSISRtsState;


/* NDF versions of the Spectrum writing */

/* Initialise the file for writing */
void acsSpecOpenTS( const char * dir,
		    unsigned int yyyymmdd,
		    unsigned int obsnum,
		    unsigned int nrecep,
		    unsigned int nsubsys,
		    const char* recepnames[],
		    int * status );

/* Write a spectrum to the file */
void acsSpecWriteTS( unsigned int subsys,
		     unsigned int nchans,
		     const float spectrum[],
		     const ACSISRtsState * state,
		     const AstFitsChan * freq,
		     int * status );

/* Close the file */
void acsSpecCloseTS( const AstFitsChan * fits[],
		   int * status );



/* INCLUDE_HDS_SPEC_H */
#endif
