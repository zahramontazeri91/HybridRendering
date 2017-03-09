# (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>
#
# LICENSE: GNU General Public License V2  
#
lib_LTLIBRARIES = libinminbnmin1.la 

libinminbnmin1_la_SOURCES= mcpoint.cc mcpoint.h markovchain.cc	\
                             markovchain.h nestedinitial.h	\
                             nestedinitial.cc nestedsampler.cc	\
                             nestedsampler.h prior_sampler.cc	\
                            prior_sampler.h

libinminbnmin1_la_CPPFLAGS = -Wall $(EIGENPKG_CFLAGS) $(BOOST_CPPFLAGS)  

libinminbnmin1_la_LDFLAGS = -version-info $(INMIN_LT_VERSION)  
