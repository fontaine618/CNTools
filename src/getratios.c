/* Jianhua Zhang. All rights reserved

*/

#include <math.h>
#include<string.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
/*parameter 1 data to be mapped to and parameter 2 segment data
   value - segment mean
   what - mean, median, max or min
   merged - vector of the same length as parameter 1 data
*/
void getratios(char **chrom1, double *start1, double *end1, int *length1, char **chrom2, double *start2,
    double  *end2, int* length2, double *value, char **what, double* merged){
	int i, j, counter, pass;
	double toreturn, medians[100], hold;
	for (i = 0; i < length1[0]; i++) {
		if(strcmp(what[0], "max") == 0){
			toreturn = -1000;
		}else if (strcmp(what[0], "min") == 0){
			toreturn = 1000;
		}else{
			toreturn = 0;
		}
		counter = 0;
		for(j = 0; j < length2[0]; j++){
		    if(strcmp(chrom1[i], chrom2[j]) == 0 && start1[i] < end2[j] && end1[i] > start2[j]){
			    if(strcmp(what[0],  "mean") == 0){
				    toreturn += value[j];
				    counter += 1;
			    }else if (strcmp(what[0], "median") == 0){
                	medians[counter] = value[j];
                    counter += 1;
				}else if(strcmp(what[0], "max") == 0){
				    if(toreturn < value[j]){
					    toreturn = value[j];
					}
			    }else if(strcmp(what[0], "min") == 0){
				    if(toreturn > value[j]){
					    toreturn = value[j];
					}
				}
			}
		}
		if(toreturn == 1000 || toreturn == -1000){
			continue;
		}else{
		    if(strcmp(what[0],  "mean") == 0 && counter != 0){
				toreturn = toreturn/counter;
		    }else if(strcmp(what[0],  "median") == 0 && counter != 0){
	            if(counter == 2){
					toreturn = (medians[0] + medians[1])/2;
				}else {
				    for(pass = 1; pass <= counter - 1; pass++){
				    	for(j = 0; j < counter - 2; j++){
				    		if(medians[j] > medians[j + 1]){
				    			hold = medians[j];
				    			medians[j] = medians[j +1];
				    			medians[j+1] = hold;
				    		}
				    	}
				    }
				    if(!(counter % 2)){
				    	toreturn = medians[counter/2];
				    }else{
				    	toreturn = medians[(counter - 1)/2];
				    }
				}
		    }
		    merged[i] = toreturn;
	    }
    }
	 return;
 }










