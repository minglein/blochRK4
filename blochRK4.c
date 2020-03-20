/**
 * Author: Mingming Wu
 * Date created: 25.11.2013
 *
 * Part Master Thesis work at ETH Zurich
 * Supervision: Markus Weiger
 *
 * Description: numnerical simulation of MR signal evolution based on Bloch differential equations
 *
 * Runge-Kutta fourth order method was implemented for robust numnerical iteration of the Bloch differential equations.
 * The increased stability is particulary advantageous when facing fast relaxation (bone signal) in a short periode of time. 
 *
 * Input: 
 *  - mxin: initial magnetization in x direction
 *  - myin: initial magnetization in y direction
 *  - mzin: initial magnetization in z direction
 *  - b1: array with magnitude of RF pulse/B1 strength 
 *  - phase: array with phase of RF pulse/B1
 *  - dt: time step between different blocks in the b1-array
 *  - froff: off resonance frequency/B0 off resonance
 *  - t1: longitudinal relaxation time
 *  - t2: transverse relaxation time
 *
 * Output: 
 *  - mx: array of magnetization in x direction during the RF pulse
 *  - my: array of magnetization in y direction during the RF pulse
 *  - mz: array of magnetization in z direction during the RF pulse
 *
 **/


#include "mex.h" 
#include <stdio.h>
#include <math.h>

#define PI 3.141592653589793

//set 3x1 vector to zero
void zerovec(double *vec)
{
*vec++=0;
*vec++=0;
*vec++=0;
}
//prepare b1
void calcBeins(double b1amp, double omega1Phase, double omega0, double *b1out)
{
*b1out++ = b1amp*cos(omega1Phase);
*b1out++ = b1amp*sin(omega1Phase);
*b1out++ = -omega0;
}
// rotation and relaxation 
void rotarela(double *Mmag, double *b1act, double t1, double t2, double *deltaMmag)
{
    if (t2>0)
    {
    deltaMmag[0] = Mmag[1]*b1act[2]-Mmag[2]*b1act[1]-Mmag[0]/t2;
    deltaMmag[1] = Mmag[2]*b1act[0]-Mmag[0]*b1act[2]-Mmag[1]/t2;
    }
    else 
    {
    deltaMmag[0] = Mmag[1]*b1act[2]-Mmag[2]*b1act[1];
    deltaMmag[1] = Mmag[2]*b1act[0]-Mmag[0]*b1act[2];
    }
            
    if (t1>0)
    {
    deltaMmag[2] = Mmag[0]*b1act[1]-Mmag[1]*b1act[0]-(Mmag[2]-1)/t1;
    }
    else
    {
    deltaMmag[2] = Mmag[0]*b1act[1]-Mmag[1]*b1act[0];    
    }
}
// runge-kutta order 4
void rkvier(double *Mmag, double *deltaMmag, double *b1act, double dt, double t1, double t2)
{
 double kx1,kx2,kx3,kx4;
 double ky1,ky2,ky3,ky4;
 double kz1,kz2,kz3,kz4;
 //double *increment;
 
 //increment = (double *) malloc(3 * sizeof(double));
 
     if (t2>0)
    {
    deltaMmag[0] = Mmag[1]*b1act[2]-Mmag[2]*b1act[1]-Mmag[0]/t2;
    deltaMmag[1] = Mmag[2]*b1act[0]-Mmag[0]*b1act[2]-Mmag[1]/t2;
    }
    else 
    {
    deltaMmag[0] = Mmag[1]*b1act[2]-Mmag[2]*b1act[1];
    deltaMmag[1] = Mmag[2]*b1act[0]-Mmag[0]*b1act[2];
    }
            
    if (t1>0)
    {
    deltaMmag[2] = Mmag[0]*b1act[1]-Mmag[1]*b1act[0]-(Mmag[2]-1)/t1;
    }
    else
    {
    deltaMmag[2] = Mmag[0]*b1act[1]-Mmag[1]*b1act[0];    
    }


 if (t2>0)
 {
     kx1=deltaMmag[0];
     kx2=(Mmag[1]+deltaMmag[1]*dt/2)*b1act[2]-(Mmag[2]+deltaMmag[2]*dt/2)*b1act[1]-(Mmag[0]+dt/2*kx1)/t2;
     kx3=(Mmag[1]+deltaMmag[1]*dt/2)*b1act[2]-(Mmag[2]+deltaMmag[2]*dt/2)*b1act[1]-(Mmag[0]+dt/2*kx2)/t2;
     kx4=(Mmag[1]+deltaMmag[1]*dt)*b1act[2]-(Mmag[2]+deltaMmag[2]*dt)*b1act[1]-(Mmag[0]+dt*kx3)/t2;
     
     ky1=deltaMmag[1];
     ky2=(Mmag[2]+deltaMmag[2]*dt/2)*b1act[0]-(Mmag[0]+deltaMmag[0]*dt/2)*b1act[2]-(Mmag[1]+dt/2*ky1)/t2;
     ky3=(Mmag[2]+deltaMmag[2]*dt/2)*b1act[0]-(Mmag[0]+deltaMmag[0]*dt/2)*b1act[2]-(Mmag[1]+dt/2*ky2)/t2;
     ky4=(Mmag[2]+deltaMmag[2]*dt)*b1act[0]-(Mmag[0]+deltaMmag[0]*dt)*b1act[2]-(Mmag[1]+dt*ky3)/t2;
  }
 else
 {
     kx1=deltaMmag[0];
     kx2=(Mmag[1]+deltaMmag[1]*dt/2)*b1act[2]-(Mmag[2]+deltaMmag[2]*dt/2)*b1act[1];
     kx3=(Mmag[1]+deltaMmag[1]*dt/2)*b1act[2]-(Mmag[2]+deltaMmag[2]*dt/2)*b1act[1];
     kx4=(Mmag[1]+deltaMmag[1]*dt)*b1act[2]-(Mmag[2]+deltaMmag[2]*dt)*b1act[1];
     
     ky1=deltaMmag[1];
     ky2=(Mmag[2]+deltaMmag[2]*dt/2)*b1act[0]-(Mmag[0]+deltaMmag[0]*dt/2)*b1act[2];
     ky3=(Mmag[2]+deltaMmag[2]*dt/2)*b1act[0]-(Mmag[0]+deltaMmag[0]*dt/2)*b1act[2];
     ky4=(Mmag[2]+deltaMmag[2]*dt)*b1act[0]-(Mmag[0]+deltaMmag[0]*dt)*b1act[2]; 
 }
 
 if (t1>0)
 {
     kz1=deltaMmag[2];
     kz2=(Mmag[0]+deltaMmag[0]*dt/2)*b1act[1]-(Mmag[1]+deltaMmag[1]*dt/2)*b1act[0]-(Mmag[2]+dt/2*kz1-1)/t1;
     kz3=(Mmag[0]+deltaMmag[0]*dt/2)*b1act[1]-(Mmag[1]+deltaMmag[1]*dt/2)*b1act[0]-(Mmag[2]+dt/2*kz2-1)/t1;
     kz4=(Mmag[0]+deltaMmag[0]*dt)*b1act[1]-(Mmag[1]+deltaMmag[1]*dt)*b1act[0]-(Mmag[2]+dt*kz3-1)/t1;
 }
 else
 {
     kz1=deltaMmag[2];
     kz2=(Mmag[0]+deltaMmag[0]*dt/2)*b1act[1]-(Mmag[1]+deltaMmag[1]*dt/2)*b1act[0];
     kz3=(Mmag[0]+deltaMmag[0]*dt/2)*b1act[1]-(Mmag[1]+deltaMmag[1]*dt/2)*b1act[0];
     kz4=(Mmag[0]+deltaMmag[0]*dt)*b1act[1]-(Mmag[1]+deltaMmag[1]*dt)*b1act[0]; 
 }
 
 deltaMmag[0]=kx1/6+kx2/3+kx3/3+kx4/6;
 deltaMmag[1]=ky1/6+ky2/3+ky3/3+ky4/6;
 deltaMmag[2]=kz1/6+kz2/3+kz3/3+kz4/6;

}

//solve bloch equation
int solvebloch(double *Mmag, double *b1, double *omega1Phase, double omega0, int ntpts, double dt, double t1, double t2, double *mx, double *my, double *mz)
{
    double deltaMmag[3];
    double b1act[3];
    double Mmagnew[3];
	double Mmagold[3];
	double beins, rfphase;
	int i;
    
	zerovec(deltaMmag);
	zerovec(b1act);
	zerovec(Mmagnew);
	
	Mmagold[0]=Mmag[0];
	Mmagold[1]=Mmag[1];
	Mmagold[2]=Mmag[2];


    for (i=0; i<ntpts; i++)
    {
    beins=*b1++;
    rfphase=*omega1Phase++;
    calcBeins(beins, rfphase, omega0, b1act);
    //rotarela(Mmagold, b1act, t1, t2, deltaMmag);
    rkvier(Mmagold, deltaMmag, b1act, dt, t1, t2);
	Mmagnew[0]=Mmagold[0]+deltaMmag[0]*dt;
	Mmagnew[1]=Mmagold[1]+deltaMmag[1]*dt;
	Mmagnew[2]=Mmagold[2]+deltaMmag[2]*dt;
    *mx=Mmagnew[0];
    *my=Mmagnew[1];
    *mz=Mmagnew[2];
	mx++;
	my++;
	mz++;
    Mmagold[0]=Mmagnew[0];
	Mmagold[1]=Mmagnew[1];
	Mmagold[2]=Mmagnew[2];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*blochRK4(mx,my,mz,omega1,phase,time,omega0off,t1,t2)*/

{
double *b1;
double *phase;    
double b0off;
double *tpts;
double dt;
double t1;
double t2;
double mxin;
double myin;
double mzin;
double *mx, *my, *mz;
double Mmag[3];

int ntpts;
int i;

mxin=*mxGetPr(prhs[0]);
myin=*mxGetPr(prhs[1]);
mzin=*mxGetPr(prhs[2]);        
b1= mxGetPr(prhs[3]);
phase = mxGetPr(prhs[4]);
dt = *mxGetPr(prhs[5]);
//off resonance frequency
b0off = *mxGetPr(prhs[6]);
//relaxation times
t1= *mxGetPr(prhs[7]);
t2= *mxGetPr(prhs[8]);

ntpts = mxGetM(prhs[3])*mxGetN(prhs[3]);



    mx = (double *)mxMalloc(ntpts*sizeof(double));
    my = (double *)mxMalloc(ntpts*sizeof(double));
    mz = (double *)mxMalloc(ntpts*sizeof(double));

/* ===== Allocate Output Magnetization vectors arrays.	*/

plhs[0] = mxCreateDoubleMatrix(ntpts,1,mxREAL);	/* Mx, output. */
plhs[1] = mxCreateDoubleMatrix(ntpts,1,mxREAL);	/* My, output. */
plhs[2] = mxCreateDoubleMatrix(ntpts,1,mxREAL);	/* Mz, output. */


mx = mxGetPr(plhs[0]);
my = mxGetPr(plhs[1]);
mz = mxGetPr(plhs[2]);


// magnetization magnitude in x, y and z direction
Mmag[0]=mxin;
Mmag[1]=myin;
Mmag[2]=mzin;

solvebloch(Mmag, b1, phase, b0off, ntpts, dt, t1, t2, mx, my, mz);
}
