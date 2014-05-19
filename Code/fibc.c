/*
 *  fibc.c
 *  
 *
 *  Created by Dongik Jang on 2011. 05. 09.
 *  Copyright 2011 Seoul National University. All rights reserved.
 *
 */

void fibc(int *n, double *outval)
{
	int i;
	double tmpval1, tmpval2[2];
	
	
	tmpval2[0] = 0;
	tmpval2[1] = 1;
	
    if(*n < 3){
		*outval = tmpval2[(*n-1)];
	}else{
		for(i=2; i<(*n); i++){
	        *outval = tmpval2[0] + tmpval2[1];
			tmpval2[0] = tmpval2[1];
	        tmpval2[1] = *outval;
		}
	}
}
