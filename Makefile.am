## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libacsisio.la
include_HEADERS = specwrite.h

libacsisio_la_SOURCES = \
	$(include_HEADERS) \
        $(C_ROUTINES)

C_ROUTINES = specwrite.c

# Testing

TESTS = spwTest
check_PROGRAMS = spwTest

spwTest_SOURCES = spwTest.c
spwTest_LDADD = libacsisio.la `ndf_link` @FCLIBS@


 
