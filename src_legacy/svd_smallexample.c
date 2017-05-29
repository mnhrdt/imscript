#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "reallysmall_svd.c"



// Transpose a square matrix in place
void trnm(double *a,int n)
{ double s,*p,*q;
  int i,j,e;
  for(i=0,e=n-1; i<n-1 ;++i,--e,a+=n+1){
    for(p=a+1,q=a+n,j=0; j<e ;++j){
      s= *p; *p++ = *q; *q=s; q+=n;
     }
   }
}

// multiply two rectangular matrices
void rmmult(double *rm,double *a,double *b,int n,int m,int l)
{ double z,*q0,*p,*q; int i,j,k;
  q0=(double *)calloc(m,sizeof(double));
  for(i=0; i<l ;++i,++rm){
    for(k=0,p=b+i; k<m ;p+=l) q0[k++]= *p;
    for(j=0,p=a,q=rm; j<n ;++j,q+=l){
      for(k=0,z=0.; k<m ;) z+= *p++ * q0[k++];
      *q=z;
     }
   }
  free(q0);
}




// example program follows



#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)



#define PMATRIX(a,m,n) print_matrix(#a,a,m,n)
static void print_matrix(char *name, double *a, int m, int n)
{
	printf("%s =", name);
	FORJ(m) {
		FORI(n)
			printf("\t%.3f", a[n*j + i]);
		printf("\n");
	}
}


// d: output size n (singular values)
// a: input size n*m (m>=n, will get trashed)
// u: output size m*m
// v: output size n*n
int old_svd(double *d, double *a, double *u, int m, double *v, int n);

// 
void svd(double *V, double *s2, double *A, int n);
int svd_regularized(double *d, double *a, double *u, int m, double *v, int n)
{
	if (n != m) exit(fprintf(stderr, "only square matrices here\n"));
	double A[n*n];
	FORI(n*n) A[i] = a[i];
	PMATRIX(A, n, n);
	svd(v, d, A, n);
	PMATRIX(A, n, n);
	PMATRIX(v, n, n);
	FORI(n) d[i] = sqrt(d[i]);
	FORJ(n) FORI(n)
		u[n*j+i] = d[j] ? A[n*j+i]/d[i] : 0;
	return 0;
}



static int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}


static void do_stuff(void)
{
	int m = 3;
	int n = 3;
	double a[n*m], d[n], u[m*m], v[n*n], dmat[n*m];

	// build and print a matrix of random numbers
	FORI(n*m) a[i] = randombounds(10, 99);
	PMATRIX(a, m, n);

	// compute the SVD
	int r = svd_regularized(d, a, u, m, v, n);
	(void)r;

	// show the results
	FORI(n*m) dmat[i] = 0;
	FORI(n) dmat[n*i+i] = d[i];
	PMATRIX(dmat, m, n);
	PMATRIX(u, m, m);
	PMATRIX(v, n, n);

	// check the reconstruction
	double ud[m*n], udv[m*n];
	rmmult(ud, u, dmat, m, m, n);
	trnm(v,n);
	rmmult(udv, ud, v, m, n, n);
	PMATRIX(udv, m, n);
}

int main(void)
{
	do_stuff();
	return 0;
}
