#include <R.h>
#include <Rmath.h>
//#include <stdlib.h>
//#include <iostream>
//#include <math.h>
using namespace std;

//#define im sqrt( complex<float>(-1) ) 
//cd /Applications/R/RScripts/ImageMoments
// R32 CMD SHLIB gpzm.cpp
//for eclipse, can use R or R64 for R64 gui
//dyn.load("gpzm.so")   in R
//out=.C("gpzmMoments",M=complex(), I=complex(), radius, theta, order, a, c=)


extern "C"  {
  
  //calculate polynomials for first order
double firstOrder(double rad,int q,int a){
	if (q >= (a/2)){
	  return(pow(rad-pow(rad,2),(a/2)) * pow(rad,q-(a/2)));
	} else {
	  return(pow(rad,q) * pow(1-rad,(a/2)));
	}	
}


  //M - moments are returned, 
  //I - image, 2-D matrix
  //radius, theta - calculated from r function, polar representation of pixel coordinates
  //order - polynomial order to calculate polynomials up to
  //a is stabilizes the result as it increases
  //optionally return c, a p by q matrix
  void gpzmMoments(double *I, double *ReM, double *ImM ,int *dim, double *radius, double *theta, int *order, double *a2, double *c, double *polyValsR, double *polyValsI, int *storePoly) {

    double a = a2[0];
    /*int *c = new (nothrow) int[order,order];
    //int *M = new (nothrow) int[order,order];
    if (c==0) {
      cout << "order is too large for memory\n";
      return;
      }*/

    //calculate constant for all p and q - p rows,q cols
    int index = 0;
    int pq=0;
    for (int q = 0; q <= order[0]; q++) {
      for (int p = q; p <= order[0]; p++) {
	pq = q*(order[0]+1)+p;
	//pochammler function
	c[pq] = sqrt((2*p+a+2)/(2*PI));	
	index = (2*q)+1;
	for(int i = 0; i < index; i++){
	   c[pq] = c[pq] * sqrt((a+1+p-q+i)/(p-q+1+i));
	}//end for i
    
      }//end for p
    }//end for q
    
    double M1=0;
    double M2=0;
    double M3=0;
    //x by y and 3 previous values needed for recurrence relation
    double *polynomial = (double*) malloc(dim[0]*dim[1]*sizeof(double));
    double *polynomial1 = (double*) malloc(dim[0]*dim[1]*sizeof(double));
    double *polynomial2 = (double*) malloc(dim[0]*dim[1]*sizeof(double));

    int xy=0;
	
    //compute all moments for all orders p and q, p rows by q cols

    //moments calculation
    for (int x = 0; x < dim[1]; x++) {
      for (int y = 0; y < dim[0]; y++) {
	xy = int(x*dim[0]+y);
	
	for (int q = 0; q <= order[0]; q++) {

	  //cmplxIm = complex <double> (cos(q*theta[xy]),-sin(q*theta[xy])) * I[xy];
	  	
	  for (int p = q; p <= order[0]; p++) {
	    pq = int(q*(order[0]+1)+p);
		
	    if(p > q+1) {
	      //constants for p and q
	      M1 = ((2*p+1+a)*(2*p+a))/((p+q+1+a)*(p-q));
	      M2 = -((p+q+1)*(a+2*p))/(p+q+a+1) + M1*((p+q)*(p-q-1))/(2*p-1+a);
	      M3 = ((p+q)*(p+q+1)*(2*p-2+a)*(2*p-1+a))/(2*(p+q+a+1)*(p+q+a)) + M2*((p+q)*(2*p-2+a))/(p+q+a) - M1*((p+q)*(p+q-1)*(p-q-2))/(2*(p+q+a));
	    }//end if

	    if (p==q) {
	      polynomial[xy] = firstOrder(radius[xy],q,a);
	      polynomial1[xy] = polynomial[xy];
    
	    } else if (p==q+1) { 
	       polynomial[xy] = ((a+3+2*q)*firstOrder(radius[xy],q+1,a)) - (2*(q+1)*firstOrder(radius[xy],q,a));
	    } else {
	      
	      //recurrence relation for polynomials
	      polynomial2[xy] = polynomial1[xy];
	      polynomial1[xy] = polynomial[xy];
	      polynomial[xy] = ((M1*radius[xy] + M2) * polynomial1[xy]) + (M3*polynomial2[xy]);
	    } //end if

		ReM[pq] = ReM[pq] + polynomial[xy] * c[pq] * cos(q*theta[xy]) * I[xy];
		ImM[pq] = ImM[pq] + polynomial[xy] * c[pq] * (-sin(q*theta[xy])) * I[xy];


	  }//end for p
	}//end for q

       	
      }//end for y
    }//end for x
    
    free(polynomial);
    free(polynomial1);
    free(polynomial2);


}//end gpzmMoments


//reconstruct image from moments
  void gpzmReconstruct(double *IR, double *ReM, double *ImM, int *d, double *radius, double *theta, int *order, double *a2, double *c) {
    int o=order[0];
    double a = a2[0];
    
  //case where order is less than 2
	
  //recalculate scaled polynomials and reconstruct image

    int xy;
    int pq;
    double M1,M2,M3;
    double *polynomial = (double*) malloc(d[0]*d[1]*sizeof(double));
    double *polynomial1 = (double*) malloc(d[0]*d[1]*sizeof(double));
    double *polynomial2 = (double*) malloc(d[0]*d[1]*sizeof(double));
    
    if (polynomial==NULL || polynomial==NULL || polynomial==NULL) {
      IR[0]=1;
      return;
    }
  
    //compute reconstructed image pixel by pixel for range of orders p and q
    for (int x = 0; x < d[1]; x++) {
      for (int y = 0; y < d[0]; y++) {
	    xy = x*d[0]+y;

	    for (int q = 0; q <= o; q++) {
	      for (int p = q; p <= o; p++) {
	       	pq = q*(o+1)+p;
	
	        if (p==q) {
		  polynomial[xy] = firstOrder(radius[xy],q,a);
		  polynomial1[xy] = polynomial[xy];
		} else if (p==q+1) { 
		  polynomial[xy] = ((a+3+2*q)*firstOrder(radius[xy],q+1,a)) - (2*(q+1)*firstOrder(radius[xy],q,a));
		} else {
		  //constants for p and q
		  M1 = ((2*p+1+a)*(2*p+a))/((p+q+1+a)*(p-q));
		  M2 = -((p+q+1)*(a+2*p))/(p+q+a+1) + M1*((p+q)*(p-q-1))/(2*p-1+a);
		  M3 = ((p+q)*(p+q+1)*(2*p-2+a)*(2*p-1+a))/(2*(p+q+a+1)*(p+q+a)) + M2*((p+q)*(2*p-2+a))/(p+q+a) - M1*((p+q)*(p+q-1)*(p-q-2))/(2*(p+q+a));
		  polynomial2[xy] = polynomial1[xy];
		  polynomial1[xy] = polynomial[xy];
		  polynomial[xy] = ((M1*radius[xy] + M2) * polynomial1[xy]) + (M3*polynomial2[xy]);
		}//end if

		if(q==0) {
		  IR[xy] = IR[xy] + ReM[pq] * c[pq] * polynomial[xy];

		} else {
		  IR[xy] = IR[xy] + 2*(ReM[pq]*cos(q*theta[xy]) + ImM[pq]*sin(q*theta[xy])) * polynomial[xy]*c[pq];
		}
	    
	  }//end for p
	}//end for q

      //	range <- max(IR)-min(IR);
      //IR <- abs(IR-min(IR))/range;
      
      }//end for y
    }//end for x
    
    free(polynomial);
    free(polynomial1);
    free(polynomial2);



}//end gpzmReconstruct


}
