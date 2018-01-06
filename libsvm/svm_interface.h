#ifndef _interface_h
#define _interface_h

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <libsvm/svm.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif


void load_model( struct svm_model * & model, char* model_file, bool predict_probability );

// Return prediction result.
double predict( svm_model* const &  model, struct svm_node* const& x,
                double feature_min[], double feature_max[], unsigned int size,
                bool predict_probability, double probability_estimates[] );



#ifdef __cplusplus
}
#endif


#endif
