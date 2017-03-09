/*
   (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

   \file inmin_mon.cc

*/

#include <stdio.h>

#include "inmin_mon.h"


void inmin_tracefn_stdout(size_t level,
                          const char *msg,
                          void *state)
{
  printf("T%i: %s\n", level, msg);
}
