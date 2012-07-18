#include <R.h>
#include <Rmath.h>


extern "C"{

void RFbyC(double *I, double *ReM, double *ImM,  double *radius, double *theta, int *D, int *Nmax, int *Mmax){
	
	// ReM=[nmax+1 mmax+1] real matrix moment
	// ImM=[nmax+1 mmax+1] imaginary matrix moment
	// radius, theta= [d1 d2] matrix
	// I=[d1 d2] matrix --- Image
	
	// Temporary values
	double Rn=0;
	
	// 
	int d1=D[0];
	int d2=D[1];
	int nmax=Nmax[0];
	int mmax=Mmax[0];
	
	int xy=0;
	int mn=0;
	
	// loop through all pixels
	for (int y=0; y<d1; y++){
	for (int x=0; x<d2; x++){
		xy=x*d1+y;
		
		for (int n=0; n<=nmax;n++){
				
				if (n==0){
					Rn=1/sqrt(radius[xy]);
				} else if ((n%2)==1){
					Rn=sqrt(2/radius[xy])*sin((n+1)*PI*radius[xy]);
				} else if ((n%2)==0){
					Rn=sqrt(2/radius[xy])*cos(n*PI*radius[xy]);
				}
				// Rn is done
				
				for (int m=0; m<=mmax;m++){
					mn=m*(nmax+1)+n;
					
					ReM[mn]= ReM[mn]+ Rn * cos(m*theta[xy]) * I[xy];
					ImM[mn]= ImM[mn]+ Rn * (-sin(m*theta[xy])) * I[xy];

				
				} // end loop of m

		} // end loop of n
	
	// end loop of pixels
	}
	}
	
}


void RFReconstructbyC(double *IR, double *ReM, double *ImM, double *radius, double *theta, int *D, int *Nmax, int *Mmax){
	// given ReM=Re(FM), ImM=Im(FM): real, and imaginary matrix of moments
	
	// given radius, theta
	// given [d1 d2]=dim(IR)=dim(radius)=dim(theta)
	
	// Temporary values
	double Rn=0;
	
	// 
	int d1=D[0];
	int d2=D[1];
	int nmax=Nmax[0];
	int mmax=Mmax[0];
	
	int xy=0;
	int mn=0;
	
	// loop through all pixels of I
	for (int y=0; y<d1; y++){
	for (int x=0; x<d2; x++){
		xy = x*d1+y;
		
		// Loop through all possible orders
		for (int n=0;n<=nmax;n++){
		
			if (n==0){
					Rn=1/sqrt(radius[xy]);
				} else if ((n%2)==1){
					Rn=sqrt(2/radius[xy])*sin((n+1)*PI*radius[xy]);
				} else if ((n%2)==0){
					Rn=sqrt(2/radius[xy])*cos(n*PI*radius[xy]);
				}
			
			// Rn is done
		
			// compute value inside summation
			
			// added<- ReMoments[n+1,1] * Rn
			IR[xy]= IR[xy]+ReM[n] * Rn;
			
			for (int m=1;m<=mmax;m++){
				mn = m*(nmax+1)+n;
				//added<- added+2*Rn*(ReMoments[n+1,m+1]*cos(m*theta)+ImMoments[n+1,m+1]*sin(m*theta))
				IR[xy]= IR[xy]+2*Rn*(ReM[mn]*cos(m*theta[xy])+ImM[mn]*sin(m*theta[xy]));
				
			}
			
		}
	}
	}
}

// close extern
}
