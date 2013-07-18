/*
 *  PES.h
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 7/16/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#ifndef __PES_H__
#define __PES_H__

namespace h2o {
    class PES {
    public:
        PES() {}
        virtual ~PES() {}
        
        virtual const char* name() const = 0;
        virtual double operator()(size_t nw, const double*) = 0; // O H H O H H
        virtual double operator()(size_t nw, const double*, double*) = 0; // O H H O H H
    };
}

#endif // __PES_H__