#include <R.h>
#include <Rmath.h>


extern "C"{

void FMbyC(double *I, double *ReM, double *ImM,  double *radius, double *theta, int *D, int *Nmax, int *Mmax){
	
	// ReM=[nmax+1 mmax+1] real matrix moment
	// ImM=[nmax+1 mmax+1] imaginary matrix moment
	// radius, theta= [d1 d2] matrix
	// I=[d1 d2] matrix --- Image
	
	// 
	int d1=D[0];
	int d2=D[1];
	int nmax=Nmax[0];
	int mmax=Mmax[0];
	
	double Rn=0;
	
	double M1=0;
    double M2=0;
    double M3=0;
    //x by y and 3 previous values needed for recurrence relation
    double *polynomial = (double*) malloc(d1*d2*sizeof(double));
    double *polynomial1 = (double*) malloc(d1*d2*sizeof(double));
    double *polynomial2 = (double*) malloc(d1*d2*sizeof(double));

    int xy=0;
	int mn=0;

	double a=0;
    //compute all moments for all orders p and q, p rows by q cols

	for (int n = 0; n <= nmax; n++) {
		
			if (n>1){
			//constants for p and q
			M1 = ((2*n+1+a)*(2*n+a))/((n+1+a)*n);
			M2 = -((n+1)*(a+2*n))/(n+a+1) + M1*(n*(n-1))/(2*n-1+a);
			M3 = (n*(n+1)*(2*n-2+a)*(2*n-1+a))/(2*(n+a+1)*(n+a)) + M2*(n*(2*n-2+a))/(n+a)-M1*(n*(n-1)*(n-2))/(2*(n+a));
		}
			//moments calculation
   		for (int x = 0; x < d2; x++) {
			for (int y = 0; y < d1; y++) {
				xy = x*d1+y;
					if (n==0) {
					polynomial[xy] = 1;//firstOrder(radius[xy],q,a);//radius^q*(1-radius)^(a/2))
					polynomial1[xy] = 1;
   
				} else if (n==1) { 
					polynomial[xy] = (a+3)*radius[xy] - 2;
				
				} else {
      
					//recurrence relation for polynomials
					polynomial2[xy] = polynomial1[xy];
					polynomial1[xy] = polynomial[xy];
					polynomial[xy] = ((M1*radius[xy] + M2) * polynomial1[xy]) + (M3*polynomial2[xy]);
	
				} //end if
				
				Rn=polynomial[xy]*(n+1)/PI;
				
				for (int m = 0; m <= mmax; m++) {
					
					mn=m*(nmax+1)+n;

					ReM[mn]= ReM[mn] + Rn * cos(m*theta[xy]) * I[xy];
					ImM[mn]= ImM[mn] + Rn * (-sin(m*theta[xy])) * I[xy];

				
				} // end for m
			}//end for y
		}//end for x
	}//end for n
   
    free(polynomial);
    free(polynomial1);
    free(polynomial2);
	
}


void FMReconstructbyC(double *IR, double *ReM, double *ImM, double *radius, double *theta, int *D, int *Nmax, int *Mmax){
	// given ReM=Re(FM), ImM=Im(FM): real, and imaginary matrix of moments
	
	// given radius, theta
	// given [d1 d2]=dim(IR)=dim(radius)=dim(theta)

	
	// 
	int d1=D[0];
	int d2=D[1];
	int nmax=Nmax[0];
	int mmax=Mmax[0];
	

	double M1=0;
    double M2=0;
    double M3=0;
    
    double a=0;

    double Rn=0;
    double Rn1=0;
    double Rn2=0;
    //double store;

    int xy=0;
	int mn=0;
	
    //compute reconstructed image pixel by pixel for range of orders n and m
    for (int x = 0; x < d2; x++) {
	for (int y = 0; y < d1; y++) {
		xy = x*d1+y;

			for (int n = 0; n <= nmax; n++) {
				/*
					//compute Rn
					store= pow(-1,n)*(n+1);
					Rn= store;
					
					if (n>=1){
						for (int s=1;s<=n;s++){
							//store= store*(-1)*radius*(n+s+1)*(n-s+1)/(s*(s+1))
							store= store*(-1)*radius[xy]*(n+s+1)*(n-s+1)/(s*(s+1));
							Rn= Rn+ store;
						}
					}
					Rn= Rn* (n+1)/PI;
				*/
				
				if (n>1){
					//constants for p and q
					M1 = ((2*n+1+a)*(2*n+a))/((n+1+a)*n);
					M2 = -((n+1)*(a+2*n))/(n+a+1) + M1*(n*(n-1))/(2*n-1+a);
					M3 = (n*(n+1)*(2*n-2+a)*(2*n-1+a))/(2*(n+a+1)*(n+a)) + M2*(n*(2*n-2+a))/(n+a)-M1*(n*(n-1)*(n-2))/(2*(n+a));
				
					Rn2 = Rn1;
					Rn1 = Rn;
					Rn = ((M1*radius[xy] + M2) * Rn1) + (M3*Rn2);
				} else if (n==0){
					Rn = 1;
					Rn1 = 1;
				} else {
					Rn = (a+3)*radius[xy] - 2;
				}
				
				IR[xy] = IR[xy] + ReM[n]* Rn;
				for (int m = 1; m <= mmax; m++) {
					mn = m*(nmax+1)+n;
					IR[xy] = IR[xy] + 2*(ReM[mn]*cos(m*theta[xy]) + ImM[mn]*sin(m*theta[xy])) * Rn;
				}
	    
			}//end for n

      
    }//end for y
    }//end for x
	
}

// close extern
}
