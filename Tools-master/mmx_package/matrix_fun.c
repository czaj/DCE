// ==================================================
// straightforward implementations of matrix multiply
// ==================================================
void multAB(double* C, double* A, double* B, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *a, *c, tmp;  
   for( i=0; i<cB; i++ ){
      for( k=0; k<cA; k++ ){
         tmp = B[i*cA+k];
         a   = A + k*rA;
         c   = C + i*rA;
         for( j=0; j<rA; j++ ){
            c[j] += tmp * a[j];
         }
      }
   }   
}
void multAtB(double* C, double* A, double* B, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *a, *b, *c;  
   for( i=0; i<cB; i++ ){
      for( k=0; k<cA; k++ ){
         a = A + k*rA;
         b = B + i*rA;
         c = C + i*cA + k;
         for( j=0; j<rA; j++ ){
            (*c) += a[j]*b[j];
         }
      }
   }
}
void multABt(double* C, double* A, double* B, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *a, *b;  
   for( j=0; j<cA; j++ ){
      a = A + j*rA;
      b = B + j*rB;      
      for( i=0; i<rB; i++ ){
         for( k=0; k<rA; k++ ){
            C[i*rA + k] += a[k]*b[i];
         }
      }
   }
}
void multAtBt(double* C, double* A, double* B, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *b, *c, tmp;  
   for( i=0; i<cA; i++ ){
      for( k=0; k<rA; k++ ){
         tmp = A[i*rA+k];
         b   = B + k*rB;
         c   = C + i;
         for( j=0; j<rB; j++ ){
            c[j*cA] += tmp * b[j];
         }
      }
   }
}

// =============================
// multiply:   C = op(A) * op(B)
// =============================
void mulMatMat(double* C, double* A, double* B,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {
#ifndef USE_BLAS // naive C implementations

   if ( (mod[0] == 'N') && (mod[1] == 'N') )
      multAB(C, A, B,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'N') )
      multAtB(C, A, B, rA, cA, cB);
   else if ( (mod[0] == 'N') && (mod[1] == 'T') )
      multABt(C, A, B, rA, cA, rB);
   else if ( (mod[0] == 'T') && (mod[1] == 'T') )
      multAtBt(C, A, B, rA, cA, rB);

#else

   // rows(Op(A)), columns(Op(A)), columns(Op(B)), rows(C)
   ptrdiff_t ropA, copA, copB, rC;  
   // can't pass consts to fortran
   ptrdiff_t rA0 = rA, rB0 = rB;    

   char modA = mod[0], modB = mod[1];
   double one = 1.0, zero = 0.0;

   if (mod[0] != 'S'){
      if ( (mod[0] == 'N') && (mod[1] == 'N') ){
         ropA  = rA;
         copA  = cA;   
         copB  = cB;
         rC    = rA;
      } else if ( (mod[0] == 'T') && (mod[1] == 'N') ){
         ropA  = cA;
         copA  = rA;   
         copB  = cB;
         rC    = cA;   
      } else if ( (mod[0] == 'N') && (mod[1] == 'T') ){
         ropA  = rA;
         copA  = cA;   
         copB  = rB;
         rC    = rA;   
      } else if ( (mod[0] == 'T') && (mod[1] == 'T') ){
         ropA  = cA;
         copA  = rA;   
         copB  = rB;
         rC    = cA;   
      }
      dgemm(&modA, &modB, &ropA, &copB, &copA, &one, A, &rA0, B, &rB0, &zero, C, &rC);
   } else {  
      char side='L', uplo = 'U';
      ropA  = rA;
      copB  = cB;     
      dsymm(&side, &uplo, &ropA, &copB,        &one, A, &rA0, B, &rB0, &zero, C, &rC);
      // why the fuck does this not work ???
   }
#endif
}

// ================================================
// square:   C = A * op(A)  or  C = 0.5*(A*B'+B*A')
// ================================================
void squareMatMat(double* C, double* A, double* B,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {
   // can't pass consts to BLAS
   ptrdiff_t rA0 = rA, cA0 = cA, rB0 = rB; 
   // rows(Op(A)), columns(Op(A)), columns(Op(B)), rows(C)
   ptrdiff_t copA, rC;  
   int i,j; 
   double temp; 

   if ( (mod[0] == 'N') ){
      copA  = cA;
      rC    = rA;
   } else {
      copA  = rA;
      rC    = cA;   
   } 

#ifndef USE_BLAS // naive C implementations

   if ((rB == 0) || (cB == 0)){  // one input  C = A*A'   
      if ( (mod[0] == 'N') )
         multABt(C, A, A, rA, cA, rA);
      else
         multAtB(C, A, A, rA, cA, cA);
   }else{
      if ( (mod[0] == 'N') )
         multABt(C, A, B, rA, cA, rB);
      else
         multAtB(C, A, B, rA, cA, cB);

      // symmetrize
      for( i=0; i<rC; i++ )
         for( j=i; j<rC; j++ ){
            temp = C[i*rC+j] + C[j*rC+i];
            C[i*rC+j] = C[j*rC+i] = 0.5*temp;   
         }
   }

#else
   char  modA = mod[0], modB = mod[1], uplo = 'U';
   double one = 1.0, zero = 0.0, half = 0.5;

   if ((!rB) && (!cB))  // one input  C = A*A'
      dsyrk(&uplo, &modA, &rC, &copA, &one, A, &rA0, &zero, C, &rC);
   else                 // two inputs C = 0.5*(A*B'+B*A')
      dsyr2k(&uplo, &modA, &rC, &copA, &half, A, &rA0, B, &rB0, &zero, C, &rC);   

   // symmetrize
   for( i=0; i<rC; i++ )
      for( j=i+1; j<rC; j++ )
          C[i*rC+j] = C[j*rC+i];

#endif
}

// =====================================
// cholesky decomposition:   C = chol(A)
// =====================================
double dot(const double* vec1, const double* vec2, const int n)
{
	int i;
	double res = 0;
	
   for( i=0; i<n; i++ )
      res += vec1[i] * vec2[i];
   
	return res;
}


int cholA(double* A, double* scratch, const int n)
{
	int i, j, rank=0;
	double tmp;

   // in-place Cholesky factorization, store 1/L(j,j) in scratch
   for( j=0; j<n; j++ )
   {
      tmp = A[j*n+j];
      if( j )
         tmp -= dot(A+j*n, A+j*n, j);

      if( tmp < 0.000000001 )
         return rank;
      else
      {
         scratch[j] = (double)(1.0/sqrt(tmp));
         rank++;
      }

      // process off-diagonal entries, modify 'A'
      for( i=j+1; i<n; i++ )
      {
         A[i*n+j] -= dot(A+i*n, A+j*n, j);
         A[i*n+j] *= scratch[j];
      }
   }

   // copy 'scratch' to diagonal of A
   for( j=0; j<n; j++ )
      A[j*n+j] = 1./scratch[j];

	return rank;
}



void chol(double* C, double* A,  const int rA) {
   int i,j, rank;
   double temp;

   // copy upper triangle into C
   for( i=0; i<rA; i++ )
      for( j=0; j<=i; j++ )
          C[i*rA+j] = A[i*rA+j];    
   
#ifndef USE_BLAS // naive C implementations
   temp = A[0];
   rank = cholA(C, A, rA);
   // chol used A as scratch, now fix it
   if (rank) A[0] = temp;
   for( i=1; i<rank; i++ )
          A[i] = A[i*rA];     
   //if decomposition failed put -1 in C
   if (rank < rA) C[0] = -1; 
#else
   ptrdiff_t rA0 = rA;
   ptrdiff_t info;
   dpotrf("U", &rA0, C, &rA0, &info );
#endif
}


// ================================
// solve linear equations   C = A\B
// ================================
void solve(double* C, double* A, double* B,
				   const int rA, const int cA, const int rB, const int cB, 
               const char *mod, double *W, const int LW, ptrdiff_t *S) {
#ifdef USE_BLAS
   int i, j, rank;
   char  uplo = 'U', side = 'L', trans = 'N', unit = 'N';
   double one = 1.0, rcond = 0.000000001;
   ptrdiff_t rA0 = rA, cA0 = cA,  cB0 = cB, Lwork=LW, info;
   ptrdiff_t rC0 = (rA>cA) ? rA : cA;
   //ptrdiff_t ptr_S = S;

   switch (mod[0]){
      case 'L':
      case 'U':
         uplo = mod[0];
         dtrsm(&side, &uplo, &trans, &unit, &rC0, &cB0, &one, A, &rA0, C, &rC0);
         break;
      case 'P':
         dposv(&uplo, &rA0, &cB0, W, &rA0, C,  &rA0, &info);// A has already been copied into W
         break;
      default:
         if (rA == cA) {
            //dgesv(&rA0, &cB0, W, &rA0, S, C,  &rA0, &info);// A has already been copied into W
            dgesv(&rA0, &cB0, W, &rA0, (ptrdiff_t*)S, C,  &rA0, &info);// A has already been copied into W
         }
         else{
            for( i=0; i<cB; i++ )
               for( j=0; j<rB; j++ )
                  C[i*rC0+j] = B[i*rB+j];
            //dgelsy(&rA0, &cA0, &cB0, A, &rA0, C, &rC0, S, &rcond, &rank, W, &Lwork, &info);
            dgelsy(&rA0, &cA0, &cB0, A, &rA0, C, &rC0,
                    (ptrdiff_t*)S, &rcond, (ptrdiff_t*)&rank, W, &Lwork, &info);
         }
   }
#endif
}
