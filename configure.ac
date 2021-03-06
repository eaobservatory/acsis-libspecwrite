dnl    Process this file with autoconf to produce a configure script
AC_REVISION($Revision$)
  
dnl    Initialisation: package name and version number
AC_INIT(acsisio, 0.7-2, timj@jach.hawaii.edu)

dnl    Require autoconf-2.50 at least
AC_PREREQ(2.50)
dnl    Require Starlink automake
AM_INIT_AUTOMAKE(1.8.2-starlink)
  
dnl    Sanity-check: name a file in the source directory -- if this
dnl    isn't found then configure will complain
AC_CONFIG_SRCDIR([specwrite.c])
 
  
dnl    Include defaults for Starlink configurations
STAR_DEFAULTS
  
dnl    Find required versions of the programs we need for configuration
AC_PROG_CC
AC_PROG_FC
AC_PROG_LIBTOOL
AC_FC_LIBRARY_LDFLAGS

dnl    Check if sys/wait.h exists and is POSIX compatible
AC_HEADER_STDC

dnl    Library functions
AC_CHECK_FUNCS([mkstemp srandomdev])

dnl    Variable sizes
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(int)

## STAR_CNF_COMPATIBLE_SYMBOLS
 
dnl    Declare the build and use dependencies for this package
STAR_DECLARE_DEPENDENCIES([build], [sae ems ndf hds cnf])
STAR_DECLARE_DEPENDENCIES([link],  [ems cnf hds ndf])
 
dnl    List the sun/ssn/... numbers which document this package and
dnl    which are present as .tex files in this directory.
dnl STAR_LATEX_DOCUMENTATION(sunxx)
 
dnl    Would like the version number as an integer
AC_DEFINE_UNQUOTED([PACKAGE_VERSION_INTEGER], $PACKAGE_VERSION_INTEGER,
                    [Integer version number, in the form major*1e6+minor*1e3+release])

dnl    If you wish to configure extra files, you can add them to this
dnl    declaration.
AC_CONFIG_FILES(Makefile component.xml)
AC_CONFIG_HEADERS( config.h )

dnl    This is the bit that does the actual work
AC_OUTPUT
