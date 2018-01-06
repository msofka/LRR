#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <libsvm/svm.h>
#include <float.h>

#include "svm_interface.h"

//:
// \file
// \brief  Interface for svmlib so that it can be called from within functions, not command prompt.
// \author Michal Sofka
// \date   March 2008



double predict( svm_model* const &  model, struct svm_node* const& x, 
                double feature_min[], double feature_max[], unsigned int size,
                bool predict_probability, double probability_estimates[] )
{
	int svm_type=svm_get_svm_type(model);
	int nr_class=svm_get_nr_class(model);
	double *prob_estimates=NULL;
	int j;

	if(predict_probability)
	{
		if (svm_type==NU_SVR || svm_type==EPSILON_SVR)
			printf("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n",svm_get_svr_probability(model));
		else
		{
			int *labels=(int *) malloc(nr_class*sizeof(int));
			svm_get_labels(model,labels);
			prob_estimates = (double *) malloc(nr_class*sizeof(double));
			//printf("labels");		
			//for(j=0;j<nr_class;j++)
			//	printf(" %d",labels[j]);
			//printf("\n");
			free(labels);
		}
	}

  double lower = -1.0;
  double upper = 1.0;


  for( unsigned int ind = 0; ind < size; ++ind ) {
    x[ind].value = lower + (upper-lower) * 
			             (x[ind].value-feature_min[ind])/
			             (feature_max[ind]-feature_min[ind]);
  }


  double v = 0.0;

	if (predict_probability && (svm_type==C_SVC || svm_type==NU_SVC))
	{
		v = svm_predict_probability(model,x,prob_estimates);
    //printf("result: %g ",v);
    for(j=0;j<nr_class;j++) {
      probability_estimates[j] = prob_estimates[j];
			//printf("%g ",prob_estimates[j]);
    }
		//printf("\n");
	}
	else
	{
		v = svm_predict(model,x);
    //printf("result: %g\n",v);
	}

	if(predict_probability)
		free(prob_estimates);

  return v;
}


void load_model( struct svm_model * & model, char* model_file, bool predict_probability )
{
  // load model
  if((model=svm_load_model( model_file ))==0)
	{
		fprintf(stderr,"can't open model file %s\n",model_file);
		exit(1);
	}
	
	if(predict_probability)
	{
		if(svm_check_probability_model(model)==0)
		{
			fprintf(stderr,"Model does not support probabiliy estimates\n");
			exit(1);
		}
	}
	else
	{
		if(svm_check_probability_model(model)!=0)
			printf("Model supports probability estimates, but disabled in prediction.\n");
	}
}
