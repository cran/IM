#include <R.h>
#include <Rmath.h>

extern "C"{

void GPZMMultiplebyC(double *I, double *ReM, double *ImM, int *dim, double *radius, double *theta, int *order, int *a2, double *c, int*NF, double *polyValsR, double *polyValsI, int *storePoly) {

    double a=double(a2[0]);
	int nF=NF[0];
	
    //calculate constant for all p and q - p rows,q cols
    int index = 0;
    int pq= 0;
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
	int pqf=0;
	int xyf=0;
	
    //compute all moments for all orders p and q, p rows by q cols
    for (int q = 0; q <= order[0]; q++) {
		for (int p = q; p <= order[0]; p++) {
			pq = q*(order[0]+1)+p;

			if(p > q+1) {
				//constants for p and q
				M1 = ((2*p+1+a)*(2*p+a))/((p+q+1+a)*(p-q));
				M2 = -((p+q+1)*(a+2*p))/(p+q+a+1) + M1*((p+q)*(p-q-1))/(2*p-1+a);
				M3 = ((p+q)*(p+q+1)*(2*p-2+a)*(2*p-1+a))/(2*(p+q+a+1)*(p+q+a)) + M2*((p+q)*(2*p-2+a))/(p+q+a)-M1*((p+q)*(p+q-1)*(p-q-2))/(2*(p+q+a));
			}

			//moments calculation
    		for (int x = 0; x < dim[1]; x++) {
				for (int y = 0; y < dim[0]; y++) {
					xy = x*dim[0]+y;

					if (p==q) {
						polynomial[xy] = pow(radius[xy],q)*pow(1-radius[xy],a/2);//firstOrder(radius[xy],q,a);//radius^q*(1-radius)^(a/2))
						polynomial1[xy] = polynomial[xy];
    
					} else if (p==(q+1)) { 
						polynomial[xy] = ((a+3+2*q)*pow(radius[xy],q+1)*pow(1-radius[xy],a/2)) - (2*(q+1)*pow(radius[xy],q)*pow(1-radius[xy],a/2));
					
					} else {
	      
						//recurrence relation for polynomials
						polynomial2[xy] = polynomial1[xy];
						polynomial1[xy] = polynomial[xy];
						polynomial[xy] = ((M1*radius[xy] + M2) * polynomial1[xy]) + (M3*polynomial2[xy]);
		
					} //end if
					
					if (storePoly[0]==1){
						index= pq*(dim[0]*dim[1])+xy;
						// store polyVals - polyVals[y][x][p][q]=[xy][pq]
						polyValsR[index]= polynomial[xy] * c[pq] * cos(q*theta[xy]);
						polyValsI[index]= polynomial[xy] * c[pq] * (-sin(q*theta[xy]));
					
						// Compute moment for each image
						for (int f=0; f<nF; f++){
							pqf=pq*nF+f;
							xyf=xy*nF+f;
							ReM[pqf]= ReM[pqf] + polyValsR[index] * I[xyf];
							ImM[pqf]= ImM[pqf] + polyValsI[index] * I[xyf];
						}
					} else {
						for (int f=0; f<nF; f++){
							pqf=pq*nF+f;
							xyf=xy*nF+f;
							ReM[pqf]= ReM[pqf] + polynomial[xy] * c[pq] * cos(q*theta[xy]) * I[xyf];
							ImM[pqf]= ImM[pqf] + polynomial[xy] * c[pq] * (-sin(q*theta[xy])) * I[xyf];
						}
					}
				}//end for y
			}//end for x

		}//end for p
    }//end for q
    
    free(polynomial);
    free(polynomial1);
    free(polynomial2);

}//end GPZMMultiplebyC

void CHFMMultiplebyC(double *I, double *ReM, double *ImM,  double *radius, double *theta, int *D, int *Nmax, int *Mmax, int *NF, double *polyValsR, double *polyValsI, int *storePoly){
	
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
	int nF=NF[0];
	
	int xy=0;
	int mn=0;
	int mnf=0;
	int xyf=0;
	int index=0;
	
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
				
				if (storePoly[0]==1){
					index= mn* d1 * d2+xy;
					polyValsR[index]= Rn * cos(m*theta[xy]);
					polyValsI[index]= Rn * (-sin(m*theta[xy]));
				
					for (int f=0; f<nF; f++){
					
						mnf=mn*nF+f;
						xyf=xy*nF+f;
						//M[n+1,m+1]<- sum(Rn*(cos(m*theta)-1i*sin(m*theta))*I)	
						ReM[mnf]=ReM[mnf] + polyValsR[index] * I[xyf];
						ImM[mnf]=ImM[mnf] + polyValsI[index] * I[xyf];
					
					}
				} else {
					for (int f=0; f<nF; f++){
					
						mnf=mn*nF+f;
						xyf=xy*nF+f;
						//M[n+1,m+1]<- sum(Rn*(cos(m*theta)-1i*sin(m*theta))*I)	
						ReM[mnf]=ReM[mnf] + Rn * cos(m*theta[xy]) * I[xyf];
						ImM[mnf]=ImM[mnf] + Rn * (-sin(m*theta[xy])) * I[xyf];
					
					}
				}
			} // end loop of m
		
		
		}// end loop of n

	
	}
	}// end loop of pixels
	
}

void RFMultiplebyC(double *I, double *ReM, double *ImM,  double *radius, double *theta, int *D, int *Nmax, int *Mmax, int *NF, double *polyValsR, double *polyValsI, int *storePoly){
	
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
	int nF=NF[0];
	
	int xy=0;
	int mn=0;
	int xyf=0;
	int mnf=0;
	int index=0;
	
	// loop through all pixels
	for (int y=0; y<d1; y++){
	for (int x=0; x<d2; x++){
		xy=x*d1+y;
		
		for (int n=0; n<=nmax;n++){
				if (radius[xy]!=0){
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
					
					if (storePoly[0]==1){
						index= mn* d1 * d2+xy;
						polyValsR[index]= Rn * cos(m*theta[xy]);
						polyValsI[index]= Rn * (-sin(m*theta[xy]));
				
						for (int f=0; f<nF; f++){
					
							mnf=mn*nF+f;
							xyf=xy*nF+f;
							//M[n+1,m+1]<- sum(Rn*(cos(m*theta)-1i*sin(m*theta))*I)	
							ReM[mnf]=ReM[mnf] + polyValsR[index] * I[xyf];
							ImM[mnf]=ImM[mnf] + polyValsI[index] * I[xyf];
					
						}
					} else {
						for (int f=0; f<nF; f++){
					
							mnf=mn*nF+f;
							xyf=xy*nF+f;
							//M[n+1,m+1]<- sum(Rn*(cos(m*theta)-1i*sin(m*theta))*I)	
							ReM[mnf]=ReM[mnf] + Rn * cos(m*theta[xy]) * I[xyf];
							ImM[mnf]=ImM[mnf] + Rn * (-sin(m*theta[xy])) * I[xyf];
					
						}
					}
				
				}// end loop of m
			}
		}
	
	
	}
	}// end loop of pixels
	
}

void FMMultiplebyC(double *I, double *ReM, double *ImM,  double *radius, double *theta, int *D, int *Nmax, int *Mmax, int *NF, double *polyValsR, double *polyValsI, int *storePoly){
	
	// ReM=[nmax+1 mmax+1] real matrix moment
	// ImM=[nmax+1 mmax+1] imaginary matrix moment
	// radius, theta= [d1 d2] matrix
	// I=[d1 d2] matrix --- Image
	
	// 
	int d1=D[0];
	int d2=D[1];
	int nmax=Nmax[0];
	int mmax=Mmax[0];
	int nF=NF[0];
	
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
	int mnf=0;
	int xyf=0;
	int index=0;
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
					polynomial[xy] = 1;
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
					
					if (storePoly[0]==1){
						index= mn* d1 * d2+xy;
						polyValsR[index]= Rn * cos(m*theta[xy]);
						polyValsI[index]= Rn * (-sin(m*theta[xy]));
				
						for (int f=0; f<nF; f++){
					
							mnf=mn*nF+f;
							xyf=xy*nF+f;
							//M[n+1,m+1]<- sum(Rn*(cos(m*theta)-1i*sin(m*theta))*I)	
							ReM[mnf]=ReM[mnf] + polyValsR[index] * I[xyf];
							ImM[mnf]=ImM[mnf] + polyValsI[index] * I[xyf];
					
						}
					} else {
						for (int f=0; f<nF; f++){
						
							mnf=mn*nF+f;
							xyf=xy*nF+f;
							//M[n+1,m+1]<- sum(Rn*(cos(m*theta)-1i*sin(m*theta))*I)	
							ReM[mnf]=ReM[mnf] + Rn * cos(m*theta[xy]) * I[xyf];
							ImM[mnf]=ImM[mnf] + Rn * (-sin(m*theta[xy])) * I[xyf];
					
						}
					}
				
				}//end for m
			}//end for y
		}//end for x

	}//end for n
    
    free(polynomial);
    free(polynomial1);
    free(polynomial2);
	
}

void GPZMreconMulti(double *IR, double *ReM, double *ImM, int *dim, double *radius, double *theta, int *Pmax, int *a1, double *c, int *NF){
	// given ReM=Re(FM), ImM=Im(FM): real, and imaginary matrix of moments
	
	// given radius, theta
	// given [d1 d2]=dim(IR)=dim(radius)=dim(theta)

	
	// 
	int d1=dim[0];
	int d2=dim[1];
	int pmax=Pmax[0];
	double a=a1[0];
	int nF= NF[0];
	
	double M1=0;
    double M2=0;
    double M3=0;

    double Rn=0;
    double Rn1=0;
    double Rn2=0;
    //double store;
	
    int xy=0;
	int pq=0;
	int xyf=0;
	int pqf=0;
	
	 //calculate constant for all p and q - p rows,q cols
    int index = 0;
    for (int q = 0; q <= pmax; q++) {
      for (int p = q; p <= pmax; p++) {
			
			pq = q*(pmax+1)+p;
			//pochammler function
			c[pq] = sqrt((2*p+a+2)/(2*PI));	
			
			index = (2*q)+1;
			for(int i = 0; i < index; i++){
				c[pq] = c[pq] * sqrt((a+1+p-q+i)/(p-q+1+i));
			}//end for i
    
      }//end for p
    }//end for q
    
    //compute reconstructed image pixel by pixel for range of orders n and m
    for (int x = 0; x < d2; x++) {
	for (int y = 0; y < d1; y++) {
		xy = x*d1+y;
		
		for (int q = 0; q<= pmax; q++){
			for (int p = q; p <= pmax; p++) {
				pq=q*(pmax+1)+p;
				
				if (p>q+1){
					//constants for p and q
					M1 = ((2*p+1+a)*(2*p+a))/((p+q+1+a)*(p-q));
					M2 = -((p+q+1)*(a+2*p))/(p+q+a+1) + M1*((p+q)*(p-q-1))/(2*p-1+a);
					M3 = ((p+q)*(p+q+1)*(2*p-2+a)*(2*p-1+a))/(2*(p+q+a+1)*(p+q+a)) + M2*((p+q)*(2*p-2+a))/(p+q+a)-M1*((p+q)*(p+q-1)*(p-q-2))/(2*(p+q+a));
			
					Rn2 = Rn1;
					Rn1 = Rn;
					Rn = (M1*radius[xy] + M2) * Rn1 + (M3*Rn2);
				} else if (p==q){
					Rn = pow(radius[xy],q)*pow(1-radius[xy],a/2);//firstOrder(radius[xy],q,a);//radius^q*(1-radius)^(a/2))
					Rn1 = Rn;
				} else {
					Rn = ((a+3+2*q)*pow(radius[xy],q+1)*pow(1-radius[xy],a/2)) - (2*(q+1)*pow(radius[xy],q)*pow(1-radius[xy],a/2));
				}
				
				for (int f=0; f<nF; f++){
					pqf=pq*nF+f;
					xyf=xy*nF+f;
					
					if (q==0){
						IR[xyf] = IR[xyf] + ReM[pqf]* c[pq] *Rn;
					} else {
						IR[xyf] = IR[xyf] + 2*(ReM[pqf]*cos(q*theta[xy]) + ImM[pqf]*sin(q*theta[xy])) * Rn *c[pq];
					}
					
				}
			}//end for p

      }// end for q
    }//end for y
    }//end for x
	
}


void RFreconMulti(double *IR, double *ReM, double *ImM, double *radius, double *theta, int *D, int *Nmax, int *Mmax, int *NF){
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
	int xyf=0;
	int mnf=0;
	
	int nF=NF[0];
	
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
			for (int f=0; f<nF; f++){
				mnf = n * nF + f;
				xyf = xy * nF + f;
				// added<- ReMoments[n+1,1] * Rn
				IR[xyf]= IR[xyf]+ReM[mnf] * Rn;
				
			}
			for (int m=1;m<=mmax;m++){
				mn = m*(nmax+1)+n;
				
				for (int f=0; f<nF; f++){
					mnf = mn * nF + f;
					xyf = xy * nF + f;
					//added<- added+2*Rn*(ReMoments[n+1,m+1]*cos(m*theta)+ImMoments[n+1,m+1]*sin(m*theta))
					IR[xyf]= IR[xyf]+2*Rn*(ReM[mnf]*cos(m*theta[xy])+ImM[mnf]*sin(m*theta[xy]));
				}
			}
			
		}
	}
	}
}

void FMreconMulti(double *IR, double *ReM, double *ImM, double *radius, double *theta, int *D, int *Nmax, int *Mmax, int *NF){
	// given ReM=Re(FM), ImM=Im(FM): real, and imaginary matrix of moments
	
	// given radius, theta
	// given [d1 d2]=dim(IR)=dim(radius)=dim(theta)

	
	// 
	int d1=D[0];
	int d2=D[1];
	int nmax=Nmax[0];
	int mmax=Mmax[0];
	
	int nF= NF[0];
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
	int xyf=0;
	int mnf=0;

    //compute reconstructed image pixel by pixel for range of orders n and m
    for (int x = 0; x < d2; x++) {
	for (int y = 0; y < d1; y++) {
		xy = x*d1+y;

			for (int n = 0; n <= nmax; n++) {
								
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
				
				for (int f=0; f<nF; f++){
					mnf = n * nF + f;
					xyf = xy * nF + f;
					IR[xyf] = IR[xyf] + ReM[mnf]* Rn;
				}
				for (int m = 1; m <= mmax; m++) {
					mn = m*(nmax+1)+n;
					for (int f=0; f<nF; f++){
						mnf = mn * nF + f;
						xyf = xy * nF + f;
						IR[xyf] = IR[xyf] + 2*(ReM[mnf]*cos(m*theta[xy]) + ImM[mnf]*sin(m*theta[xy])) * Rn;
					}
				}
	    
			}//end for n

      
    }//end for y
    }//end for x
	
}

void CHFreconMulti(double *IR, double *ReM, double *ImM, double *radius, double *theta, int *D, int *Nmax, int *Mmax, int *NF){
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
	int nF=NF[0];
	
	int xy=0;
	int mn=0;
	int xyf=0;
	int mnf=0;
	
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
		
			for (int f=0; f<nF; f++){
					mnf = n * nF + f;
					xyf = xy * nF + f;
					IR[xyf] = IR[xyf] + ReM[mnf]* Rn;
			}
			for (int m = 1; m <= mmax; m++) {
				mn = m*(nmax+1)+n;
				for (int f=0; f<nF; f++){
					mnf = mn * nF + f;
					xyf = xy * nF + f;
					IR[xyf] = IR[xyf] + 2*(ReM[mnf]*cos(m*theta[xy]) + ImM[mnf]*sin(m*theta[xy])) * Rn;
				}
			}

		}
	}
	}
}


}// end extern
