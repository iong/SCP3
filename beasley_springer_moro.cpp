/*
 *  beasley_springer_moro.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/18/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <cmath>

#include "beasley_springer_moro.h"

void beasley_springer_moro(int nb_dims, double u[], double x[])
{
    int i, j;
    double r;
    
    double a[] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
    double b[] = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
    double c[] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
        0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
        0.0000321767881768, 0.0000002888167364, 0.0000003960315187};
    
    for (j=0; j < nb_dims; j++) {
        double yj = u[j] - 0.5;
        
        if (fabs(yj) < 0.42) {
            r = yj*yj;
            x[j] = yj*(((a[3]*r + a[2])*r + a[1])*r + a[0])/((((b[3]*r + b[2])*r + b[1])*r + b[0])*r + 1);
        }
        else {
            if (yj > 0) {
                r = log(-log(1-u[j]));
            }
            else if (yj < 0) {
                r = log(-log(u[j]));
            }
            x[j] = c[0] + r*(c[1] + r*(c[2] + r*(c[3] + r*(c[4] + r*(c[5] + r*(c[6] + r*(c[7] + r*c[8])))))));
            
            if (yj < 0) {
                x[j] = -x[j];
            }
        }
    }
}
