/************************************************************************************************/
/* SVD Using Jacobis Rotations									*/
/*												*/
/* Compile: g++ -O3 SVD.cpp -o SVD								*/
/* Arguments:											*/
/*												*/
/*	M = # of columns									*/
/*	N = # of Rows										*/
/*												*/
/*	Matrix must be squared (M=N)								*/
/*												*/
/*	-t = print out Timing and # of Iterations						*/
/*	-p = print out Results (U, S, V)							*/
/*	-d0 = Generate the Octave files for debug and verify correctness			*/
/*  -dm = Generate the Matlab files for debug and verify correctness				*/
/*												*/
/* Use:	./SVD M N -t -p -d									*/
/*												*/
/* All arguments aren't important, just M and N. If you want, is possible to do 		*/
/* ./SVD M N -t and only print out the timing. As well you can use ./SVD M N -d for debug.	*/
/************************************************************************************************/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>
#include <omp.h>
#define epsilon 1.e-8
#define MIN(x, y) ((x) < (y) ? (x) : (y))

using namespace std;

template <typename T> double sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}

int main (int argc, char* argv[]){
  int M,N;
	bool octave; 
  string T,P,Db;
  M = atoi(argv[1]);
  N = atoi(argv[2]);

  double elapsedTime,elapsedTime2;
  timeval start,end,end2;

  if(argc < 4){
	  cout<<"Please input the size of Matrix and at least one of the options: -t -p -do/-dm";
	  return 0;
  }

 
  if(M != N){
	  cout<<"Error: Matrix must be square";
	  return 0;
  }
  
  if(argc > 3){

    T = argv[3];
    if(argc > 4){
      P = argv[4];
      if(argc > 5){
        Db = argv[5];
		if(Db=="-do")
		{
			cout << "Using Octave" << endl;
			octave = true;
		}
		else if(Db=="-dm")
		{
			cout << "Using Matlab" << endl;
			octave = false;
		}
		else{
			cout << "Invalid debug option %s" << Db <<  endl;
			return 1;
		}
      }
    }
  }
 // cout<<T<<P<<endl;
  
  double **U,**V, *S,**U_t, **V_t, **A;
  double alpha, beta, gamma, c, zeta, t,s,sub_zeta, converge;

  int acum = 0;
  int temp1, temp2;
  converge = 1.0;

  U = new double*[N];
  V = new double*[N];
  U_t = new double*[N];
  V_t = new double*[N];
  A = new double*[N];
  S = new double[N];

  for(int i =0; i<N; i++){
	U[i] = new double[N];
 	V[i] = new double[N];
	U_t[i] = new double[N];
	V_t[i] = new double[N];
	A[i] = new double[N];
  }


  //Read from file matrix, if not available, app quit
  //Already transposed

  ifstream matrixfile("matrix");
  if(!(matrixfile.is_open())){
    cout<<"Error: file not found"<<endl;
    return 0;
  }

  for(int i = 0; i < M; i++){
    for(int j =0; j < N; j++){

      matrixfile >> U_t[i][j];
    }
  }

  matrixfile.close();

 
  for(int i=0; i<M;i++){
    for(int j=0; j<N;j++){

      if(i==j){
        V_t[i][j] = 1.0;
      }
      else{
        V_t[i][j] = 0.0;
      }
    }
  }


    //Store A for debug purpouse
  

   for(int i=0; i<M;i++){
      for(int j=0; j<N;j++){

       A[i][j] = U_t[j][i];
      }
    }


  

  
  /* SVD using Jacobi algorithm (Sequencial)*/

   gettimeofday(&start, NULL);

   double conv;
   while(converge > epsilon){ 		//convergence
    converge = 0.0;	
   		
    acum++;				//counter of loops

    for(int i = 1; i<M; i++){
//#pragma omp parallel private(zata, t, c, s) shared(i)
	int j;

	//for(int k = 0; k < N; k++){
	//alpha = alpha + (U_t[i][k] * U_t[i][k]);
//	}
//#pragma omp parallel for num_threads(4) private(j, alpha, beta, gamma, zeta, t, c, s) reduction(max:converge)
      for(j = 0; j<i; j++){
          alpha = 0.0;
          beta = 0.0;
          gamma = 0.0;
	int k, kk;
//	double temp[3] = {0.0, 0.0, 0.0};
	//#pragma omp parallel num_threads(4) for private(kk, k, alpha, beta, gamma) shared(N)
//	#pragma omp parallel for private(k, kk) num_threads(4) shared(N) schedule(static, N/4) 
          for(k = 0; k<N ; k+=1){
//		for(k = kk; k < MIN(N, kk + 16); k++){ 
            alpha = alpha + (U_t[i][k] * U_t[i][k]);
            beta = beta + (U_t[j][k] * U_t[j][k]);
            gamma = gamma + (U_t[i][k] * U_t[j][k]);


//		temp[0] += (U_t[i][k] * U_t[i][k]);
//		temp[1] += (U_t[j][k] * U_t[j][k]);
//		temp[2] += (U_t[i][k] * U_t[j][k]);
          }
//	}
//  }	
//	alpha = temp[0];
//	beta = temp[1];
//	gamma = temp[2];
          converge = max(converge, abs(gamma)/sqrt(alpha*beta));	//compute convergence
	  								//basicaly is the angle
									//between column i and 
          zeta = (beta - alpha) / (2.0 * gamma);
          t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));        //compute tan of angle
          c = 1.0 / (sqrt (1.0 + (t*t)));				//extract cos
          s = c*t;							//extrac sin
 

	  //Apply rotations on U and V
//	int k;
//	#pragma omp parallel for private(t, k) num_threads(4) shared(N) schedule(static, N/4)
  	  for(k=0; k<N; k+=1){
//		for(k = kk; k < MIN(N, kk + 16); k++){
	            t = U_t[i][k];
        	    U_t[i][k] = c*t - s*U_t[j][k];
       		    U_t[j][k] = s*t + c*U_t[j][k];

	            t = V_t[i][k];
         	    V_t[i][k] = c*t - s*V_t[j][k];
        	    V_t[j][k] = s*t + c*V_t[j][k];

          }
//	}
      }
    }
 }


  //Create matrix S
int i, ii, j;
#pragma omp parallel for private(i, ii, j, t) num_threads(4) schedule(static, M/4)
//#pragma omp parallel for
  for(ii =0; ii<M; ii+=16){
	for( i = ii; i < MIN(M, ii+16); i++){
    	t=0;
	//int j,jj;

//#pragma omp parallel for private(j) num_threads(4) shared(N) schedule(static, N/4) reduction(+:t) 
    for(j=0; j<N;j+=1){
//	for(int j = jj; j < MIN(N, jj + 16); j++){
 	     t=t + pow(U_t[i][j],2);
  //  }
}
    t = sqrt(t);
//#pragma omp parallel for
	#pragma omp parallel for num_threads(4) schedule(static, N/4)
    for(j=0; j<N;j+=1){
//	for(int j = jj; j < MIN(N, jj + 16); j++){
   	   U_t[i][j] = U_t[i][j] / t;
     	 if(i == j){
     	   S[i] = t;
    	  }
	}
    //  }
    }
}
  gettimeofday(&end, NULL);
 /************************************************************/

 /* Develop SVD Using OpenMP */



// fix final result

  for(int i =0; i<M; i++){
    
    for(int j =0; j<N; j++){

      U[i][j] = U_t[j][i];
      V[i][j] = V_t[j][i];
      
    }
    
  }


  //Output time and iterations


  if(T=="-t" || P =="-t"){
    cout<<"iterations: "<<acum<<endl;
    elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
    cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;


  }


  // Output the matrixes for debug
  if(T== "-p" || P == "-p"){
  cout<<"U"<<endl<<endl;
  for(int i =0; i<M; i++){
    for(int j =0; j<N; j++){

      cout<<U[i][j]<<"  ";
    }
    cout<<endl;
  }

  cout<<endl<<"V"<<endl<<endl;
  for(int i =0; i<M; i++){
    for(int j =0; j<N; j++){

      cout<<V[i][j]<<"  ";
    }
    cout<<endl;
  }

  cout<<endl<<"S"<<endl<<endl;
  for(int i =0; i<M; i++){
    for(int j =0; j<N; j++){

       if(i==j){  cout<<S[i]<<"  ";}
	
       else{
	       cout<<"0.0  ";
       }
    }
    cout<<endl;
  }

  }

  //Generate Octave files for debug purpouse
   if(Db == "-do" || T == "-do" || P == "-do" || Db == "-dm" || T == "-dm" || P == "-dm"){


    ofstream Af;
    //file for Matrix A
	if(octave){
		Af.open("matrixA.mat"); 
		Af<<"# Created from debug\n# name: A\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";
	}else{
		Af.open("matrixA.dat");
	}

    for(int i = 0; i<M;i++){
      for(int j =0; j<N;j++){
        Af<<" "<<A[i][j];
      }
      Af<<"\n";
    }
    
    Af.close();

    ofstream Uf;

    //File for Matrix U
	if(octave)
	{
		Uf.open("matrixUcpu.mat");
		Uf<<"# Created from debug\n# name: Ucpu\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";
	}else{
		Uf.open("matrixUcpu.dat");
	}
    for(int i = 0; i<M;i++){
      for(int j =0; j<N;j++){
        Uf<<" "<<U[i][j];
      }
      Uf<<"\n";
    }
    Uf.close();

    ofstream Vf;
    //File for Matrix V
	if(octave){
	    Vf.open("matrixVcpu.mat");
		Vf<<"# Created from debug\n# name: Vcpu\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";
	}else{
		Vf.open("matrixVcpu.dat");
	}
    for(int i = 0; i<M;i++){
      for(int j =0; j<N;j++){
        Vf<<" "<<V[i][j];
      }
      Vf<<"\n";
    }
    

    Vf.close();

    ofstream Sf;
    //File for Matrix S
    if(octave){
		Sf.open("matrixScpu.mat");
		Sf<<"# Created from debug\n# name: Scpu\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";
	}else{
		Sf.open("matrixScpu.dat");
	}
    
    for(int i = 0; i<M;i++){
      for(int j =0; j<N;j++){
        if(i == j){
         Sf<<" "<<S[i];

        }

        else{
          Sf<<" 0.0";
        }
      }
      Sf<<"\n";
    }
    

    Sf.close();


 }

   delete [] S;
   for(int i = 0; i<N;i++){
	   delete [] A[i];
	   delete [] U[i];
	   delete [] V[i];
	   delete [] U_t[i];
	   delete [] V_t[i];
   }

  return 0;
}






