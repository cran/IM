#include <R.h>
#include <Rmath.h>


extern "C"{

void CHFMbyC(double *I, double *ReM, double *ImM,  double *radius, double *theta, int *D, int *Nmax, int *Mmax){
	
	// ReM=[nmax+1 mmax+1] real matrix moment
	// ImM=[nmax+1 mmax+1] imaginary matrix moment
	// radius, theta= [d1 d2] matrix
	// I=[d1 d2] matrix --- Image
	
	// Temporary values
	double consT=0;
	double Rn=0;
	double temP=0;
	
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
		xy = x*d1+y;
		
		consT= sqrt(8/PI)*pow(-1+1/radius[xy],1/4);
		
		// Loop through all possible orders
		for (int n=0;n<=nmax;n++){
			
			// compute Rn(r)
			
			// n=0,n=1; compute Rn=store=const*(2*(2*radius-1))^n
			temP=consT*pow(2*(2*radius[xy]-1),n);
			Rn=temP;
			
			if (n>=2){
				// Compute Rn with n>=2
				for (int k=1;k<=(n/2);k++){
					// store=store*(-1)*(n-2*(k-1))*(n-2*(k-1)-1)/(n-(k-1))/((k-1)+1)*pow((2*(2*radius-1)),(-2))
					// Rn[i][j]=Rn[i][j]+store
					temP= temP*(-1)*(n-2*(k-1))*(n-2*(k-1)-1)/(n-(k-1))/((k-1)+1)*pow(2*(2*radius[xy]-1),-2);
					Rn= Rn+temP;
				}
			}
			// Rn is done computed
			
			for (int m=0; m<=mmax;m++){
				
				mn = m*(nmax+1)+n;

				ReM[mn]=ReM[mn]+ Rn * cos(m*theta[xy]) * I[xy];
				ImM[mn]=ImM[mn]+ Rn * (-sin(m*theta[xy])) * I[xy];

			}
		
		// end loop of n
		}

	// end loop of pixels
	}
	}
	
}


void CHFMReconstructbyC(double *IR, double *ReM, double *ImM, double *radius, double *theta, int *D, int *Nmax, int *Mmax){
	// given ReM=Re(CHFM), ImM=Im(CHFM): real, and imaginary matrix of moments
	
	// given radius, theta
	// given [d1 d2]=dim(IR)=dim(radius)=dim(theta)
	
	// Temporary values
	double consT=0;
	double Rn=0;
	double temP=0;
	
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
		consT= sqrt(8/PI)*pow(-1+1/radius[xy],1/4);
		
		// Loop through all possible orders
		for (int n=0;n<=nmax;n++){
		
			// n=0,n=1; compute Rn=store=const*(2*(2*radius-1))^n
			temP=consT*pow(2*(2*radius[xy]-1),n);
			Rn=temP;
			
			if (n>=2){
				// Compute Rn with n>=2
				for (int k=1;k<=(n/2);k++){
					// store=store*(-1)*(n-2*(k-1))*(n-2*(k-1)-1)/(n-(k-1))/((k-1)+1)*pow((2*(2*radius-1)),(-2))
					// Rn[i][j]=Rn[i][j]+store
					temP= temP*(-1)*(n-2*(k-1))*(n-2*(k-1)-1)/(n-(k-1))/((k-1)+1)*pow(2*(2*radius[xy]-1),-2);
					Rn= Rn+temP;
				}
			}
		
			// compute value inside summation
			// added<- ReMoments[n+1,1] * Rn
			IR[xy]= IR[xy]+ ReM[n] * Rn;
			
			for (int m=1;m<=mmax;m++){
			
				mn= m*(nmax+1)+n;
				//added<- added+2*Rn*(ReMoments[n+1,m+1]*cos(m*theta)+ImMoments[n+1,m+1]*sin(m*theta))
				IR[xy]= IR[xy]+2*Rn*(ReM[mn]*cos(m*theta[xy])+ImM[mn]*sin(m*theta[xy]));
				
			}
			
		}
	}
	}
}

// close extern
}
