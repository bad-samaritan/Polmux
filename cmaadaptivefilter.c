#include "mex.h"

/* Author: Massimiliano Salsi, 2009			*/
/* University of Parma, Italy				*/

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

/* Multiply and Accumulate a vector */
double vmac(double *a, double *b, int n)
{
    double z = 0;
    int i=0;
    for(i=0;i<n;i++)
        z += a[i] * b[i];
    return z;
}

double errorfun(double a, double b, double M)
{
    /* errorfun = (M - a*a - b*b) */;
    return (M - a*a - b*b);
}

/* Updating of coefficients */
void updatecoeff(double *hr,double *hi, int Ntap,
                 double *xr,double *xi, int Ndim,
                 double outr, double outi, double mu, double R)
{
    double k = mu * errorfun(outr,outi,R);
    int i;
    for(i=0;i<Ntap;i++)
    {
        *(hr+i) += k * (outr * xr[i] + outi * xi[i]);
        *(hi+i) += k * (outi * xr[i] - outr * xi[i]);
        *(hr+i+Ntap) += k * (outr * xr[i+Ndim] + outi * xi[i+Ndim]);
        *(hi+i+Ntap) += k * (outi * xr[i+Ndim] - outr * xi[i+Ndim]); 
    }
}

/* Core CMA adaptive filter function */
void cmafilter(double *xr, double *xi, int Ndim, 
               double *h1r,double *h1i,double *h2r,double *h2i, int Ntap,
               double mu,double *R,
               double *yr,double *yi,
               bool dontskip)
{
    int i;
    int k = ((Ntap - 1) / 2) % 2;
    int dimY = Ndim-Ntap+1;
    for(i=0;i<dimY;i++) /* dimY */
    {
        /* Calculating output field */
        
        /* Output first column real part yr = Real( x * h1) */
        *(yr+i) = vmac(xr+i,h1r,Ntap) - vmac(xi+i,h1i,Ntap) /* 1st filter column */
               +vmac(xr+i+Ndim,h1r+Ntap,Ntap)-vmac(xi+i+Ndim,h1i+Ntap,Ntap); /* 2nd filter column */;
        /* Output first column imag part yi = Imag( x * h1) */
        *(yi+i) = vmac(xi+i,h1r,Ntap) + vmac(xr+i,h1i,Ntap) /* 1st filter column */
               +vmac(xi+i+Ndim,h1r+Ntap,Ntap)+vmac(xr+i+Ndim,h1i+Ntap,Ntap); /* 2nd filter column */;
        /* Output second column real part yr = Real( x * h2) */
        *(yr+i+dimY) = vmac(xr+i,h2r,Ntap) - vmac(xi+i,h2i,Ntap) /* 1st filter column */
               +vmac(xr+i+Ndim,h2r+Ntap,Ntap)-vmac(xi+i+Ndim,h2i+Ntap,Ntap); /* 2nd filter column */;
        /* Output second column imag part yi = Imag( x * h1) */
        *(yi+i+dimY) = vmac(xi+i,h2r,Ntap) + vmac(xr+i,h2i,Ntap) /* 1st filter column */
               +vmac(xi+i+Ndim,h2r+Ntap,Ntap)+vmac(xr+i+Ndim,h2i+Ntap,Ntap); /* 2nd filter column */;
        
        /* Updating filter coefficients */
        
        if (  dontskip || ( i % 2 == k) )
        {
            updatecoeff(h1r,h1i,Ntap,xr+i,xi+i,Ndim,*(yr+i),*(yi+i),mu,*R);
            updatecoeff(h2r,h2i,Ntap,xr+i,xi+i,Ndim,*(yr+i+dimY),*(yi+i+dimY),mu,*(R+1));
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray *prhs[])
{
    double *xr  = mxGetPr(prhs[0]);     /* Input field vector, real part */
    double *xi  = mxGetPi(prhs[0]);     /* Input field vector, imag part */
    double *h1r = mxGetPr(prhs[1]);     /* 1st output filter, real part  */
    double *h1i = mxGetPi(prhs[1]);     /* 1st output filter, imag part  */
    double *h2r = mxGetPr(prhs[2]);     /* 2nd output filter, real part  */
    double *h2i = mxGetPi(prhs[2]);     /* 2nd output filter, imag part  */
    double Ntap = mxGetScalar(prhs[3]); /* Length of filters, real scalar*/
    double mu   = mxGetScalar(prhs[4]); /* Algo constant, real scalar    */
    double *R   = mxGetPr(prhs[5]);     /* Convergence radii, real vector*/
	double sps  = mxGetScalar(prhs[6]); /* Samples x symbol, 1 or 2      */
    
    int Mdim     = mxGetM(prhs[0]);     /* Length of input field vector  */
    int Npol     = mxGetN(prhs[0]);     /* Should be 2                   */
    int Mfilter1 = mxGetM(prhs[1]);     /* Should be == Ntap             */
    int Nfilter1 = mxGetN(prhs[1]);     /* Should be 2                   */
    int Mfilter2 = mxGetM(prhs[2]);     /* Should be == Ntap             */
    int Nfilter2 = mxGetN(prhs[2]);     /* Should be 2                   */
    
    bool dontskip;
    
    /* Will be better to add some checks and error:
     * Ntap == Mfilter1 == Mfilter2 == ODD INTEGER
     * Npol == Nfilter1 == Nfilter2 == 2 */
    if ((int)Ntap % 2 == 0)
        mexErrMsgTxt("Ntaps should be an ODD INTEGER.");
    
    /* Chechking the value of samples per symbol    */
    if ((int)sps == 1)
        dontskip = 1;
    else
    {
        if ((int)sps == 2)
        {
            dontskip = 0;
        }
        else
        {
            mexErrMsgTxt("Samples x symbol should be either 1 or 2.");
        }
    }
        
    /* Chech if the vector that are supposed to be complex are really complex */
    /* If X, H1 or H2 are purely real numbers (which is often the case for H1 */
    /* and H2) then matlab does not allocates the memory for the imaginary    */
    /* part and the program will segfault if it tries to access it            */
    
    if (xi==NULL)
    {   /* If necessary allocates memory for the imaginary part of X */
        xi = mxCalloc(Mdim*Npol, sizeof(double));
        mxSetPi( (mxArray *)prhs[0], xi);
    }
    if (h1i==NULL)
    {   /* If necessary allocates memory for the imaginary part of H1 */
        h1i = mxCalloc(Ntap*Npol,sizeof(double));
        mxSetPi( (mxArray *)prhs[1], h1i);
    }
    if (h2i==NULL)
    {   /* If necessary allocates memory for the imaginary part of H2 */
        h2i = mxCalloc(Ntap*Npol,sizeof(double));
        mxSetPi( (mxArray *)prhs[2], h2i);
    }
        
    
	/* Next line allocates memory for a (Mdim-Ntap+1) x Npol complex vector */
    plhs[0] = mxCreateDoubleMatrix(Mdim-Ntap+1,Npol,mxCOMPLEX);
    
	double *yr = mxGetPr(plhs[0]);     /* Output field vector, real part */
    double *yi = mxGetPi(plhs[0]);     /* Output field vector, imag part */
    
    cmafilter(xr,xi,Mdim,h1r,h1i,h2r,h2i,Ntap,mu,R,yr,yi,dontskip);
   
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 
    double *y = mxGetPr(plhs[1]);
    *y = 0;
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *z = mxGetPr(plhs[2]);
    *z = 0;

    return;
}
