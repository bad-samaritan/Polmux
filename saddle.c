/***************************************************************************
 *            saddle.c
 *
 * This file is the .c version of saddle.m. See saddle.m for help.
 * Compile this file with 
 *
 *	mex saddle.c
 *
 * If the mex file is created successfully, matlab discards saddle.m and
 * works with the mex implementation that is much more faster.
 *
 *  Author Marco Bertolini, 2009
 *  Email bertolini@tlc.unipr.it
 ****************************************************************************/

/*
 *    This file is part of Optilux, the optical simulator toolbox.
 *    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
 *			 
 *    Optilux is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 * 
 *    Optilux is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#define SADTOL 1e-5  /*tolerance of the Newton's method*/
#define SADACC 1e-3  /*accuracy if MAXNX is reached*/
#define MAXNX  100   /*maximum number of iterations*/
#define TOLDEN 1e-9  /*avoid division by very small numbers*/

/*Evaluate the saddlepoint u0, see [1]. nsw = 0: ok. nsw = 1: no
  saddlepoint. nsw = 2: inaccurate saddlepoint.
*/

int sign (double x)
{
	if (x < 0)  {return -1;};
	if (x == 0) {return 0;};
	if (x > 0)  {return 1;};
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *sgn  = mxGetPr(prhs[0]); 
	double *xi   = mxGetPr(prhs[1]);
	double *dx   = mxGetPr(prhs[2]);
	int     nbit = mxGetM(prhs[0]);
	double  dsg2 = mxGetScalar(prhs[3]);
	double  qsg2 = mxGetScalar(prhs[4]);
	double  qsg4 = mxGetScalar(prhs[5]);
	double  vars = mxGetScalar(prhs[6]);
	double *ld   = mxGetPr(prhs[7]);
	int     len_ld = mxGetM(prhs[7]);
	double *ldex = mxGetPr(prhs[8]);
	double *ld2  = mxGetPr(prhs[9]);
	double *b2   = mxGetPr(prhs[10]);
	double  pol2 = mxGetScalar(prhs[11]);
	 
	
	plhs[0] = mxCreateDoubleMatrix(1,nbit,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,nbit,mxREAL);
	
	int *nsw = (int*)mxGetPr(plhs[1]);        /*1: no saddle point. 2: inaccurate saddle point.*/
	double ubl=0;
    double ubr=0;
    double *u0 = mxGetPr(plhs[0]);
    int k;
    
    /* Get the value of eps */
    double eps=mxGetEps();
    
    for( k=0; k<nbit; k++ )
    {
        nsw[k] = 0;
    /*Set up the initial guess*/
        u0[k] = (dx[k]+sgn[k]*sqrt(dx[k]*dx[k]+4*vars))/(2*vars);  /*(B.11) of [1]*/
        
        if (u0[k] > 0)
        {
            double ubi = 1/(dsg2*ldex[1]);
            ubl = 0;
            ubr = ubi;
            if (u0[k] >= ubi) { u0[k] = 0.99*ubi;}  /*see comment after (B.11) in [1]*/
        }
        else
        {
            ubl = -1E300;
            ubr = 0;
            if (ldex[0] < 0)
            {
                ubl = 1/(dsg2*ldex[0]);
                if (u0[k] <= ubl) {u0[k] = 0.99*ubl;} /*see comment after (B.11) in [1]*/
            }
        }
        
        
    /*Search the saddle point with the Newton's method*/
        int    nx = 0;
        double du = 1E300;
        double u1,u2 = 0;
        
        while ((fabs(du) >= SADTOL) && (nx <= MAXNX))
        {
            
            int i;
            double umbs,fden,fnum,sum1=0,sum2=0;
            
            double ds2u0 = dsg2*u0[k];
            double d1 = -(xi[k]+1/u0[k]);
            double d2 = 1/(u0[k]*u0[k]);
            for(i=0;i<len_ld;i++)
            {
                umbs=1-ds2u0*ld[i];
                fden=umbs*umbs;
                fnum=b2[i+len_ld*k]/ld[i] + dsg2*ld[i]*umbs*pol2;
                sum1+=fnum/fden;
                fden=fden*umbs;
                fnum=qsg2*b2[i+len_ld*k]+qsg4*pol2*ld2[i]*umbs;
                sum2+=fnum/fden;
            }
            
            d1 = d1 + sum1;  /*first derivative, (B.8) of [1]*/
            d2 = d2 + sum2;  /* second derivative, (B.10) of [1]*/
            
            if (d2 < eps)
            {
                if (nsw[k] == 0)
                {
                    nsw[k] = 1;
                    return;
                }
            }
            
            nx = nx+1;
            u2 = u1;    /*saddle point at cycle nx-2*/
            u1 = u0[k];    /*saddle point at cycle nx-1*/
            du = -d1/d2;
            u0[k] = u1+du; /*saddle point */
            
            while ((u0[k] <= ubl) || (u0[k] >= ubr))  /* Out of bounds? */
            {
                du = 0.5*du;		    /* correct. */
                u0[k] = u1 + du;
            }
            
        }
        
        if (nx >= MAXNX && (sign(u2-u1) != sign(u1-u0[k])))
        {
            if (fabs((u1-u0[k])/u0[k]) >= SADACC && nsw[k] == 0)
            {
                nsw[k] = 2;
            }
        }
        double den = u0[k]-2*u1+u2;
        if ((nx >= 2) && (den > TOLDEN))
        {
            du = u0[k]-u1;
            u0[k] = u0[k]-du*du/den;	/* Aitken's acceleration */
        }
    }
    
}
