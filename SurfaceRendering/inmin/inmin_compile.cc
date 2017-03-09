/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
   (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

*/

#include <iostream>

#include "inmin_compile.h"

#include "libtcc.h"

// Note that definition of inmin_s1d_data must be identical to the
// definition in the main compiler (see inmin_compile.h)
char inmin_preamble[] = \
"double  acos(double);"
"double asin(double);"
"double atan(double);"
"double atan2(double);"
"double ceil(double);"
"double cos(double);"
"double cosh(double);"
"double exp(double);"
"double fabs(double);"
"double floor(double);"
"double fmod(double);"
"double frexp(double);"
"double ldexp(double);"
"double log(double);"
"double log10(double);"
"double modf(double);"
"double pow(double);"
"double sin(double);"
"double sinh(double);"
"double sqrt(double);"
"double tan(double);"
"double tanh(double);"
"typedef struct {\n" 
"unsigned M;\n"
"const double *xa;\n"
"const double *ya;\n"
  "} inmin_s1d_data;\n";

typedef struct {
    unsigned M;
    const double *xa;
    const double *ya;
} inminfit_data;


static void inmin_tcc_reporterr(void *d,
                                    const char *msg)
{
    std::cerr<<msg<<std::endl;
}

inmin_f_res inmin_s1d_fn_c(const char *f)
{
    std::string prog=std::string(inmin_preamble)+std::string("void inminfit_ff(const double *p,  double *res,       inmin_s1d_data *data){\n")+
        std::string("for(int i=0; i<data->M; ++i){")+
        std::string("const double x=data->xa[i];\n")+
        std::string("const double y=("+std::string(f)+");\n")+
        std::string("res[i]=data->ya[i]-y;}}");
    
    TCCState *s=tcc_new();
    tcc_enable_debug(s);
    
    
    tcc_set_output_type(s, 
                        TCC_OUTPUT_MEMORY);


    tcc_set_error_func(s, 
                       NULL,
                       inmin_tcc_reporterr);
    
    tcc_compile_string(s, 
                       prog.c_str());

    
    // The function code will be pointed by theis pointer
    inmin_f_res ff_func;
    
    if(tcc_relocate(s)==-1)
        {
            std::cerr<<"Relocations failed!";
        };
    
    // Do not delete s!

    ff_func = (inmin_f_res)( tcc_get_symbol(s, "inminfit_ff"));
    return ff_func;
}

