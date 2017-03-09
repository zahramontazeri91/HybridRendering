# (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>
#
# LICENSE: GNU General Public License V2  
#
SUBDIRS = bnmin1

lib_LTLIBRARIES = libinmin.la 

libinmin_la_SOURCES= inmin_lm.cc inmin_lm.h \
                     inmin_mon.cc inmin_mon.h	\
                     eigen_lm.cc eigen_utils.hpp  \
		     eigen_qo.cc eigen_qo.hpp \
		     inmin_compile.h \
		     inmin_nested.cc inmin_nested.h

libinmin_la_CFLAGS = -Wall 
libinmin_la_CPPFLAGS = -Wall $(EIGENPKG_CFLAGS) $(BOOST_CPPFLAGS) 
libinmin_la_LIBADD = $(TINYCC_LDFLAGS) bnmin1/libinminbnmin1.la

libinmin_la_LDFLAGS = -version-info $(INMIN_LT_VERSION) 

if TINYCC_OPT
   libinmin_la_SOURCES += inmin_compile.cc 
   libinmin_la_CPPFLAGS += $(TINYCC_CPPFLAGS)
   libinmin_la_LDFLAGS += -ltcc -ldl
else
   libinmin_la_SOURCES += inmin_compile_stubs.cc 
endif

include_HEADERS = inmin_lm.h inmin_mon.h inmin_compile.h inmin_nested.h
