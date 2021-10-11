/*     
 * File:   cel_vector.h -- see headers for details and function's comments
 * Project CRESTA (see details on https://cresta-project.eu) 
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#include "cel_vector.h"
#include <stdlib.h>
#include <stdio.h>
cel_vector* cel_vector_new(INDEX vector_num, int output_on)
{
	INDEX ii;
	void* temp_alloc;
	int err_code;
	cel_vector*  vector;
	err_code=0;
	vector = Nullptr(cel_vector);     
	err_code += posix_memalign(&temp_alloc, (size_t)POSIX_MEM_ALIGN, sizeof(cel_vector)*(vector_num));
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on>1)
			fprintf(stderr,"ERROR in cel_vector_new -> not enought memory\n");
		return Nullptr(cel_vector);
	}
	vector =  (cel_vector*)temp_alloc;
	for(ii=0; ii<vector_num;ii++)
	{
		vector[ii].length = 0;
		vector[ii].values = Nullptr(REAL);
	}
	return vector;
}

int cel_vector_allocate(INDEX elements_num, cel_vector* vector, int output_on)
{
	int err_code;
	void* temp_alloc;

	err_code=0;

	err_code += posix_memalign(&temp_alloc, (size_t)64, sizeof(REAL)*(elements_num));
	vector->values = temp_alloc;
	
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on==1)
			fprintf(stderr,"ERROR in cel_vector_allocate -> not enought memory\n");
		return err_code;
	}
	vector->length = elements_num;

	return err_code;
}

void cel_vector_free(cel_vector* vector, INDEX vector_num)
{
	INDEX ii;
	for(ii=0; ii<vector_num;ii++)
	{
		free(vector[ii].values);
		vector[ii].values = Nullptr(REAL);
		vector[ii].length=0;
	}
}

void cel_vector_print(cel_vector* vector, int output_on)
{
	INDEX ii;
	if(output_on)
	{
		fprintf(stderr,"VECTOR (length="INDEX_FORMAT"):\n", vector->length);
		
		for(ii=0; ii<vector->length; ii++)				
		{
			fprintf(stderr,"\n("INDEX_FORMAT"): "REAL_FORMAT";\n",ii,vector->values[ii]);
		}

		fprintf(stderr,"\nVECTOR END\n");
	}
}

