// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define SOR_EXTRAPOLATION_PARAMETER 1.9
#define INPUT_PRESMOOTHING_SIGMA 0.8



// this function is like "malloc", but it returns always a valid pointer
static void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p)
		exit(fprintf(stderr, "out of memory\n"));
	return p;
}

#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

/**
 *
 * Compute the gradient of an image using centered differences
 *
 * For computing the derivatives near the boundary of the image, it is assumed
 * that the pixels outside the image domain take the value of the nearest pixel
 * inside the image domain.
 *
 */
static void compute_gradient_using_centered_differences(
	const float *input, // input image
	float *dx,          // computed x derivative
	float *dy,          // computed y derivative
	const int nx,       // image width
	const int ny        // image height
)
{
	// compute gradient in the central body of the image
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i = 1; i < ny-1; i++)
	{
		for(int j = 1; j < nx-1; j++)
		{
			const int k = i * nx + j;
			dx[k] = 0.5*(input[k+1] - input[k-1]);
			dy[k] = 0.5*(input[k+nx] - input[k-nx]);
		}
	}

	// compute gradient in the first and last rows
	for(int j = 1; j < nx-1; j++)
	{
		dx[j] = 0.5*(input[j+1] - input[j-1]);
		dy[j] = 0.5*(input[j+nx] - input[j]);

		const int k = (ny - 1) * nx + j;

		dx[k] = 0.5*(input[k+1] - input[k-1]);
		dy[k] = 0.5*(input[k] - input[k-nx]);
	}

	// compute gradient in the first and last columns
	for(int i = 1; i < ny-1; i++)
	{
		const int p = i * nx;
		dx[p] = 0.5*(input[p+1] - input[p]);
		dy[p] = 0.5*(input[p+nx] - input[p-nx]);

		const int k = (i+1) * nx - 1;

		dx[k] = 0.5*(input[k] - input[k-1]);
		dy[k] = 0.5*(input[k+nx] - input[k-nx]);
	}

	// compute the gradient in the corners
	dx[0] = 0.5*(input[1] - input[0]);
	dy[0] = 0.5*(input[nx] - input[0]);

	dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
	dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

	dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
	dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

	dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
	dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);

}


/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
static void gaussian_smoothing_in_place(
	float *I,             // input/output image
	const int xdim,       // image width
	const int ydim,       // image height
	const double sigma    // Gaussian sigma
)
{
	const int boundary_condition = DEFAULT_BOUNDARY_CONDITION;
	const int window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

	const double den  = 2*sigma*sigma;
	const int   size = (int) (window_size * sigma) + 1 ;
	const int   bdx  = xdim + size;
	const int   bdy  = ydim + size;

	if (boundary_condition && size > xdim) {
		fprintf(stderr, "GaussianSmooth: sigma too large (s=%g,w=%d) {%d}\n", sigma, size, xdim);
		abort();
	}

	// compute the coefficients of the 1D convolution kernel
	double B[size];
	for(int i = 0; i < size; i++)
		B[i] = 1 / (sigma * sqrt(2.0 * M_PI)) * exp(-i * i / den);

	// normalize the 1D convolution kernel
	double norm = 0;
	for(int i = 0; i < size; i++)
		norm += B[i];
	norm *= 2;
	norm -= B[0];
	for(int i = 0; i < size; i++)
		B[i] /= norm;

	// convolution of each line of the input image
	double *R = xmalloc((size + xdim + size)*sizeof(double));

	for (int k = 0; k < ydim; k++)
	{
		int i, j;
		for (i = size; i < bdx; i++)
			R[i] = I[k * xdim + i - size];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for(i = 0, j = bdx; i < size; i++, j++)
				R[i] = R[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + size-i];
				R[j] = I[k * xdim + xdim-i-1];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + xdim-size+i];
				R[j] = I[k * xdim + i];
			}
			break;
		}

		for (i = size; i < bdx; i++)
		{
			double sum = B[0] * R[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * ( R[i-j] + R[i+j] );
			I[k * xdim + i - size] = sum;
		}
	}

	// convolution of each column of the input image
	double *T = xmalloc((size + ydim + size)*sizeof(double));

	for (int k = 0; k < xdim; k++)
	{
		int i, j;
		for (i = size; i < bdy; i++)
			T[i] = I[(i - size) * xdim + k];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for (i = 0, j = bdy; i < size; i++, j++)
				T[i] = T[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for (i = 0, j = bdy; i < size; i++, j++) {
				T[i] = I[(size-i) * xdim + k];
				T[j] = I[(ydim-i-1) * xdim + k];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for( i = 0, j = bdx; i < size; i++, j++) {
				T[i] = I[(ydim-size+i) * xdim + k];
				T[j] = I[i * xdim + k];
			}
			break;
		}

		for (i = size; i < bdy; i++)
		{
			double sum = B[0] * T[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * (T[i-j] + T[i+j]);
			I[(i - size) * xdim + k] = sum;
		}
	}

	free(R);
	free(T);
}



#define BICUBIC_BOUNDARY_NEUMANN 0
#define BICUBIC_BOUNDARY_PERIODIC 1
#define BICUBIC_BOUNDARY_SYMMETRIC 2

#define DEFAULT_BICUBIC_BOUNDARY_CONDITION BICUBIC_BOUNDARY_NEUMANN

/**
  *
  * Neumann boundary condition test
  *
**/
static int neumann_bc(int x, int nx, bool *out)
{
	if (x < 0)
	{
	    x = 0;
	    *out = true;
	}
	else if (x >= nx)
	{
	    x = nx - 1;
	    *out = true;
	}

	return x;
}

/**
  *
  * Periodic boundary condition test
  *
**/
static int periodic_bc(int x, int nx, bool *out)
{
	if (x < 0)
	{
		const int n   = 1 - (int)(x/(nx+1));
		const int ixx = x + n * nx;

		x =   ixx% nx;
		*out = true;
	}
	else if (x >= nx)
	{
		x = x % nx;
		*out = true;
	}

	return x;
}


/**
  *
  * Symmetric boundary condition test
  *
**/
static int symmetric_bc(int x, int nx, bool *out)
{
	if (x < 0)
	{
		const int top = nx - 1;
		const int xx = -x;
		const int n  = (int)(xx/top) % 2;

		if ( n ) x = top - ( xx % top );
		else x = xx % top;
		*out = true;
	}
	else if ( x >= nx )
	{
		const int top = nx - 1;
		const int n = (int)(x/top) % 2;

		if ( n ) x = top - ( x % top );
		else x = x % top;
		*out = true;
	}

	return x;
}


/**
  *
  * Cubic interpolation in one dimension
  *
  * This function evaluates a polynomial P(x) which has
  * the property that P(i) = v[i+1] for i=-1,0,1,2
  *
  * The coefficients of such a polynomial can be easily calculated by hand.
  *
**/
static double cubic_interpolation_cell (
	double v[4],  // interpolation values
	double x      // point to be interpolated
)
{
	return  v[1] + 0.5 * x * (v[2] - v[0] +
		x * (2.0 *  v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] +
		x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}


/**
  *
  * Bicubic interpolation in two dimensions
  *
  * This is a separable polynomial P(x,y)=f(x)*g(y) 
  * that interpolates the 16 points of a 4x4 cell.
  *
**/
static double bicubic_interpolation_cell (
	double p[4][4], // array containing the interpolation values
	double x,       // x position to be interpolated
	double y        // y position to be interpolated
)
{
	double v[4];
	v[0] = cubic_interpolation_cell(p[0], y);
	v[1] = cubic_interpolation_cell(p[1], y);
	v[2] = cubic_interpolation_cell(p[2], y);
	v[3] = cubic_interpolation_cell(p[3], y);
	return cubic_interpolation_cell(v, x);
}

/**
  *
  * Compute the bicubic interpolation of a point in an image.
  * Detect if the point goes outside the image domain.
  *
**/
static float bicubic_interpolation_at(
	const float *input, //image to be interpolated
	const float  uu,    //x component of the vector field
	const float  vv,    //y component of the vector field
	const int    nx,    //image width
	const int    ny,    //image height
	bool         border_out //if true, return zero outside the region
)
{
	const int boundary_condition = DEFAULT_BICUBIC_BOUNDARY_CONDITION;
	const int sx = (uu < 0) ? -1: 1;
	const int sy = (vv < 0) ? -1: 1;

	int x, y, mx, my, dx, dy, ddx, ddy;
	bool out[1] = {false};

	//apply the corresponding boundary conditions
	switch(boundary_condition)
	{
	case BICUBIC_BOUNDARY_NEUMANN:
		x   = neumann_bc((int) uu, nx, out);
		y   = neumann_bc((int) vv, ny, out);
		mx  = neumann_bc((int) uu - sx, nx, out);
		my  = neumann_bc((int) vv - sx, ny, out);
		dx  = neumann_bc((int) uu + sx, nx, out);
		dy  = neumann_bc((int) vv + sy, ny, out);
		ddx = neumann_bc((int) uu + 2*sx, nx, out);
		ddy = neumann_bc((int) vv + 2*sy, ny, out);
		break;

	case BICUBIC_BOUNDARY_PERIODIC:
		x   = periodic_bc((int) uu, nx, out);
		y   = periodic_bc((int) vv, ny, out);
		mx  = periodic_bc((int) uu - sx, nx, out);
		my  = periodic_bc((int) vv - sx, ny, out);
		dx  = periodic_bc((int) uu + sx, nx, out);
		dy  = periodic_bc((int) vv + sy, ny, out);
		ddx = periodic_bc((int) uu + 2*sx, nx, out);
		ddy = periodic_bc((int) vv + 2*sy, ny, out);
		break;

	case BICUBIC_BOUNDARY_SYMMETRIC:
		x   = symmetric_bc((int) uu, nx, out);
		y   = symmetric_bc((int) vv, ny, out);
		mx  = symmetric_bc((int) uu - sx, nx, out);
		my  = symmetric_bc((int) vv - sx, ny, out);
		dx  = symmetric_bc((int) uu + sx, nx, out);
		dy  = symmetric_bc((int) vv + sy, ny, out);
		ddx = symmetric_bc((int) uu + 2*sx, nx, out);
		ddy = symmetric_bc((int) vv + 2*sy, ny, out);
		break;
	}

	if(*out && border_out)
		return 0.0;

	else
	{
		//obtain the interpolation points of the image
		const float p11 = input[mx  + nx * my];
		const float p12 = input[x   + nx * my];
		const float p13 = input[dx  + nx * my];
		const float p14 = input[ddx + nx * my];

		const float p21 = input[mx  + nx * y];
		const float p22 = input[x   + nx * y];
		const float p23 = input[dx  + nx * y];
		const float p24 = input[ddx + nx * y];

		const float p31 = input[mx  + nx * dy];
		const float p32 = input[x   + nx * dy];
		const float p33 = input[dx  + nx * dy];
		const float p34 = input[ddx + nx * dy];

		const float p41 = input[mx  + nx * ddy];
		const float p42 = input[x   + nx * ddy];
		const float p43 = input[dx  + nx * ddy];
		const float p44 = input[ddx + nx * ddy];

		//create array
		double pol[4][4] = {
			{p11, p21, p31, p41},
			{p12, p22, p32, p42},
			{p13, p23, p33, p43},
			{p14, p24, p34, p44}
		};

		//return interpolation
		return bicubic_interpolation_cell(pol, uu-x, vv-y);
	}
}


/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
static void bicubic_interpolation_warp(
	const float *input,  //image to be warped
	const float *u,      //x component of the vector field
	const float *v,      //y component of the vector field
	float       *output, //warped output image with bicubic interpolation
	const int    nx,     //image width
	const int    ny,     //image height
	bool         border_out//if true, put zeros outside the region
)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i = 0; i < ny; i++)
		for(int j = 0; j < nx; j++)
		{
			const int   p  = i * nx + j;
			const float uu = (float) (j + u[p]);
			const float vv = (float) (i + v[p]);

			//obtain the bicubic interpolation at position (uu, vv)
			output[p] = bicubic_interpolation_at(input,
					uu, vv, nx, ny, border_out);
		}
}




#define ZOOM_SIGMA_ZERO 0.6

/**
  *
  * Compute the size of a zoomed image from the zoom factor
  *
**/
static void zoom_size(
	int nx,             // width of the orignal image
	int ny,             // height of the orignal image
	int *nxx,           // width of the zoomed image
	int *nyy,           // height of the zoomed image
	float factor        // zoom factor between 0 and 1
)
{
	*nxx = round(nx * factor);
	*nyy = round(ny * factor);
}

/**
  *
  * Downsample an image
  *
**/
static void zoom_out(
	const float *I,          // input image
	float *Iout,             // output image
	const int nx,            // image width
	const int ny,            // image height
	const float factor       // zoom factor between 0 and 1
)
{
	// temporary working image
	float *Is = xmalloc(nx * ny * sizeof(float));
	for(int i = 0; i < nx * ny; i++)
		Is[i] = I[i];

	// compute the size of the zoomed image
	int nxx, nyy;
	zoom_size(nx, ny, &nxx, &nyy, factor);

	// compute the Gaussian sigma for smoothing
	const float sigma = ZOOM_SIGMA_ZERO * sqrt(1.0/(factor*factor) - 1.0);

	// pre-smooth the image
	gaussian_smoothing_in_place(Is, nx, ny, sigma);

	// re-sample the image using bicubic interpolation
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i1 = 0; i1 < nyy; i1++)
	for (int j1 = 0; j1 < nxx; j1++)
	{
		const float i2  = (float) i1 / factor;
		const float j2  = (float) j1 / factor;

		float g = bicubic_interpolation_at(Is, j2, i2, nx, ny, false);
		Iout[i1 * nxx + j1] = g;
	}

	free(Is);
}


/**
  *
  * Function to upsample the image
  *
**/
static void zoom_in(
	const float *I, // input image
	float *Iout,    // output image
	int nx,         // width of the original image
	int ny,         // height of the original image
	int nxx,        // width of the zoomed image
	int nyy         // height of the zoomed image
)
{
	// compute the zoom factor
	const float factorx = ((float)nxx / nx);
	const float factory = ((float)nyy / ny);

	// re-sample the image using bicubic interpolation
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i1 = 0; i1 < nyy; i1++)
	for (int j1 = 0; j1 < nxx; j1++)
	{
		float i2 =  (float) i1 / factory;
		float j2 =  (float) j1 / factorx;

		float g = bicubic_interpolation_at(I, j2, i2, nx, ny, false);
		Iout[i1 * nxx + j1] = g;
	}
}



/**
 *
 *  Function to compute the SOR iteration at a given position
 *  (SOR = Successive Over-Relaxation)
 *
 *  This function corresponds to lines 8-10 on the pseudocode
 *  written on the associated article for the procedure
 *  "horn_schunck_optical_flow"
 *
 */
static float sor_iteration(
	const float *Au, // constant part of the numerator of u
	const float *Av, // constant part of the numerator of v
	const float *Du, // denominator of u
	const float *Dv, // denominator of v
	const float *D,  // constant part of the numerator
	float 	    *u,  // x component of the flow
	float 	    *v,  // y component of the flow
	const float  al, // alpha smoothness parameter
	const int    p,  // current position
	const int    p1, // up-left neighbor
	const int    p2, // up-right neighbor
	const int    p3, // bottom-left neighbor
	const int    p4, // bottom-right neighbor
	const int    p5, // up neighbor
	const int    p6, // left neighbor
	const int    p7, // bottom neighbor
	const int    p8  // right neighbor
	// scheme of the neighbor positions:
	//
	// p1 p5 p2
	// p6 p  p8
	// p3 p7 p4
)
{
	// set the SOR extrapolation parameter
	const float w = SOR_EXTRAPOLATION_PARAMETER;

	// compute the divergence
	const float ula = 1./12. * (u[p1] + u[p2] + u[p3] + u[p4]) +
			  1./6.  * (u[p5] + u[p6] + u[p7] + u[p8]);
	const float vla = 1./12. * (v[p1] + v[p2] + v[p3] + v[p4]) +
			  1./6.  * (v[p5] + v[p6] + v[p7] + v[p8]);

	// store the previous values
	const float uk = u[p];
	const float vk = v[p];

	// update the flow
	u[p] = (1.0 - w) * uk + w * (Au[p] - D[p] * v[p] + al * ula) / Du[p];
	v[p] = (1.0 - w) * vk + w * (Av[p] - D[p] * u[p] + al * vla) / Dv[p];

	// return the convergence error
	return (u[p] - uk) * (u[p] - uk) + (v[p] - vk) * (v[p] - vk);
}



/**
 *
 *  Horn & Schunck method for optical flow estimation at a single scale
 *
 */
static void horn_schunck_optical_flow(
	const float *I1,             // source image
	const float *I2,             // target image
	float       *u,              // x component of optical flow
	float       *v,              // y component of optical flow
	const int    nx,             // image width
	const int    ny,             // image height
	const float  alpha,          // smoothing parameter
	const int    nwarps,         // number of warpings per scale
	const float  epsilon,        // stopping criterion threshold
	const int    nmaxiter,       // maximum number of iterations
	const bool   verbose         // switch on messages
)
{
	if (verbose) fprintf(stderr, "Single-scale Horn-Schunck of a %dx%d "
		       "image\n\ta=%g nw=%d eps=%g mi=%d v=%d\n", nx, ny,
			alpha, nwarps, epsilon, nmaxiter, verbose);

	const int   npixels = nx * ny;
	const float alpha2  = alpha * alpha;

	//allocate memory
	int sf = sizeof(float);
	float *I2x  = xmalloc(npixels * sf); // x derivative of I2
	float *I2y  = xmalloc(npixels * sf); // y derivative of I2
	float *I2w  = xmalloc(npixels * sf); // warping of I2
	float *I2wx = xmalloc(npixels * sf); // warping of I2x
	float *I2wy = xmalloc(npixels * sf); // warping of I2y
	float *Au   = xmalloc(npixels * sf); // constant part of numerator of u
	float *Av   = xmalloc(npixels * sf); // constant part of numerator of v
	float *Du   = xmalloc(npixels * sf); // denominator of u
	float *Dv   = xmalloc(npixels * sf); // denominator of v
	float *D    = xmalloc(npixels * sf); // common numerator of u and v

	// compute the gradient of the second image
	compute_gradient_using_centered_differences(I2, I2x, I2y, nx, ny);

	// iterative approximation to the Taylor expansions
	for (int n = 0; n < nwarps; n++)
	{
		if(verbose) fprintf(stderr, "Warping %d:", n);

		// warp the second image and its derivatives
		bicubic_interpolation_warp(I2,  u, v, I2w,  nx, ny, true);
		bicubic_interpolation_warp(I2x, u, v, I2wx, nx, ny, true);
		bicubic_interpolation_warp(I2y, u, v, I2wy, nx, ny, true);

		// store the constant parts of the system
		// (including the starting values of (u,v) at this warp,
		// denoted by u^n and v^n in the pseudocode of the article)
		for(int i = 0; i < npixels; i++)
		{
			const float I2wl = I2wx[i] * u[i] + I2wy[i] * v[i];
			const float dif  = I1[i] - I2w[i] + I2wl;

			Au[i] = dif * I2wx[i];
			Av[i] = dif * I2wy[i];
			Du[i] = I2wx[i] * I2wx[i] + alpha2;
			Dv[i] = I2wy[i] * I2wy[i] + alpha2;
			D[i]  = I2wx[i] * I2wy[i];
		}

		// counter for the loop below (named "r" on article pseudocode)
		int niter = 0;

		// this variable starts with a dummy value to enter the loop
		float stopping_criterion = epsilon + 1;

		// iterations of the SOR numerical scheme
		while (stopping_criterion > epsilon && niter < nmaxiter)
		{
			niter++;
			stopping_criterion = 0;

			//process the central part of the optical flow
#ifdef _OPENMP
			#pragma omp parallel for reduction(+:stopping_criterion)
#endif
			for(int i = 1; i < ny-1; i++)
			for(int j = 1; j < nx-1; j++)
			{
				const int k = i * nx + j;
				stopping_criterion += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx-1, k-nx+1, k+nx-1,
						k+nx+1, k-nx, k-1, k+nx, k+1
						);
			}

			// process the first and last rows
			for(int j = 1; j < nx-1; j++)
			{
				// first row
				int k = j;
				stopping_criterion += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-1, k+1, k+nx-1, k+nx+1,
						k, k-1, k+nx, k+1
						);

				// last row
				k = (ny-1) * nx + j;
				stopping_criterion += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx-1, k-nx+1, k-1, k+1,
						k-nx, k-1, k, k+1
						);
			}

			// process the first and last columns
			for(int i = 1; i < ny-1; i++)
			{
				// first column
				int k = i * nx;
				stopping_criterion += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx, k-nx+1, k+nx, k+nx+1,
						k-nx, k, k+nx, k+1
						);

				// last column
				k = (i+1) * nx - 1;
				stopping_criterion += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx-1, k-nx, k+nx-1, k+nx,
						k-nx, k-1, k+nx, k
						);
			}

			// process the corners
			// up-left corner
			stopping_criterion += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					0, 0, 1, nx, nx+1,
					0, 0, nx, 1
					);

			// up-right corner
			int k = nx - 1;
			stopping_criterion += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					k, k-1, k, k+nx-1, k+nx,
					k, k-1, k+nx, k
					);

			// bottom-left corner
			k = (ny-1) * nx;
			stopping_criterion += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					k, k-nx, k-nx+1,k, k+1,
					k-nx, k, k, k+1
					);

			// bottom-right corner
			k = ny * nx - 1;
			stopping_criterion += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					k, k-1, k, k-nx-1, k-nx,
					k-nx, k-1, k, k
					);

			stopping_criterion = sqrt(stopping_criterion / npixels);
		}

		if(verbose)
			fprintf(stderr, "Iterations %d (%g)\n",
				       	niter, stopping_criterion);
	}

	// free the allocated memory
	free(I2x);
	free(I2y);
	free(I2w);
	free(I2wx);
	free(I2wy);
	free(Au);
	free(Av);
	free(Du);
	free(Dv);
	free(D);
}

// compute the largest number of an array
static float max_element(const float *x, int n)
{
	int r = 0;
	for (int i = 1; i < n; i++)
		if (x[i] > x[r])
			r = i;
	return x[r];
}

// compute the smallest number of an array
static float min_element(const float *x, int n)
{
	int r = 0;
	for (int i = 1; i < n; i++)
		if (x[i] < x[r])
			r = i;
	return x[r];
}


/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
static void image_normalization(
	const float *I1,      // first image
	const float *I2,      // second image
	float       *I1n,     // first normalized image
	float       *I2n,     // second normalized image
	int          npixels  // size of each image
)
{
	// find the max and min of both images
	const float max1 = max_element(I1, npixels);
	const float max2 = max_element(I2, npixels);
	const float min1 = min_element(I1, npixels);
	const float min2 = min_element(I2, npixels);

	// obtain the absolute max and min
	const float max = max1 > max2 ? max1 : max2;
	const float min = min1 < min2 ? min1 : min2;
	const float den = max - min;

	if(den > 0)
		// normalize both images
		for(int i = 0; i < npixels; i++)
		{
			I1n[i] = 255.0 * (I1[i] - min) / den;
			I2n[i] = 255.0 * (I2[i] - min) / den;
		}

	else
		// copy the original images
		for(int i = 0; i < npixels; i++)
		{
			I1n[i] = I1[i];
			I2n[i] = I2[i];
		}
}


/**
 *
 *  Procedure to handle the pyramidal approach.
 *  This procedure relies on the previous functions to calculate
 *  large optical flow fields using a pyramidal scheme.
 *
 */
 void horn_schunck_pyramidal(
	const float *I1,              // source image
	const float *I2,              // target image
	float       *u,               // x component of optical flow
	float       *v,               // y component of optical flow
	const int    nx,              // image width
	const int    ny,              // image height
	const float  alpha,           // smoothing weight
	      int    nscales,         // number of scales
	const float  zfactor,         // zoom factor
	const int    nwarps,          // number of warpings per scale
	const float  epsilon,         // stopping criterion threshold
	const int    nmaxiter,         // maximum number of iterations
	const bool   verbose          // switch on messages
)
{
	float N = 1 + log(hypot(nx, ny)/16)/log(1/zfactor);
	if (N < nscales)
		nscales = N;

	if (verbose) fprintf(stderr, "Multiscale Horn-Schunck of a %dx%d pair"
			"\n\ta=%g ns=%d zf=%g nw=%d eps=%g mi=%d\n", nx, ny,
			alpha, nscales, zfactor, nwarps, epsilon, nmaxiter);

	int npixels = nx * ny;

	float *I1s[nscales];
	float *I2s[nscales];
	float *us[nscales];
	float *vs[nscales];
	int nxx[nscales];
	int nyy[nscales];


	I1s[0] = xmalloc(npixels * sizeof(float));
	I2s[0] = xmalloc(npixels * sizeof(float));

	// normalize the finest scale images between 0 and 255
	image_normalization(I1, I2, I1s[0], I2s[0], npixels);

	// presmoothing the finest scale images
	gaussian_smoothing_in_place(I1s[0], nx, ny, INPUT_PRESMOOTHING_SIGMA);
	gaussian_smoothing_in_place(I2s[0], nx, ny, INPUT_PRESMOOTHING_SIGMA);

	us[0] = u;
	vs[0] = v;
	nxx[0] = nx;
	nyy[0] = ny;

	// create the scales
	for(int s = 1; s < nscales; s++)
	{
		zoom_size(nxx[s-1], nyy[s-1], nxx+s, nyy+s, zfactor);

		const int npixels_s = nxx[s] * nyy[s];
		I1s[s] = xmalloc(npixels_s * sizeof(float));
		I2s[s] = xmalloc(npixels_s * sizeof(float));
		us[s] =  xmalloc(npixels_s * sizeof(float));
		vs[s] =  xmalloc(npixels_s * sizeof(float));

		// compute the zoom from the previous finer scale
		zoom_out(I1s[s-1], I1s[s], nxx[s-1], nyy[s-1], zfactor);
		zoom_out(I2s[s-1], I2s[s], nxx[s-1], nyy[s-1], zfactor);
	}

	// initialize the flow
	for (int i = 0; i < nxx[nscales-1] * nyy[nscales-1]; i++)
	{
		us[nscales-1][i] = 0;
		vs[nscales-1][i] = 0;
	}

	// pyramidal approximation to the optic flow
	for(int s = nscales-1; s >= 0; s--)
	{
		if(verbose)
			fprintf(stderr, "Scale: %d %dx%d\n", s, nxx[s], nyy[s]);

		// compute the optical flow at this scale
		horn_schunck_optical_flow(
			I1s[s], I2s[s], us[s], vs[s], nxx[s], nyy[s],
			alpha, nwarps, epsilon, nmaxiter, verbose
		);

		// if this was the last scale, finish now
		if (!s) break;

		// otherwise, upsample the optical flow

		// zoom the optic flow for the next finer scale
		zoom_in(us[s], us[s-1], nxx[s], nyy[s], nxx[s-1], nyy[s-1]);
		zoom_in(vs[s], vs[s-1], nxx[s], nyy[s], nxx[s-1], nyy[s-1]);

		// scale the optic flow with the appropriate zoom factor
		for(int i = 0; i < nxx[s-1] * nyy[s-1]; i++)
		{
			us[s-1][i] *= 1.0 / zfactor;
			vs[s-1][i] *= 1.0 / zfactor;
		}
	}

	// free the allocated memory
	free(I1s[0]);
	free(I2s[0]);
	for(int i = 1; i < nscales; i++)
	{
		free(I1s[i]);
		free(I2s[i]);
		free(us[i]);
		free(vs[i]);
	}
}
