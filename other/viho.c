// Homography viewer.
//
// This program loads a color image and allows the user to drag the four
// corners of the image around the window, thus deforming the image by an
// homography.
//
//
// Usage:
//
// 	viho image.png
//
//
// Mouse controls:
//
// 	1. Drag the control points with the left mouse button to change their
// 	position in the window
//
//	2. Drag the control points with the right mouse button to change their
//	position in the image
//
//	3. Mouse wheel for zoom-in and zoom-out
//
//
// Keys:
//
//	q	exit the viewer
//	c	reset viewer and center image
//
// 	+	zoom-in
// 	-	zoom-out
// 	ARROWS	move the view
//
// 	0	use nearest-neighbor interpolation
// 	1	use linear interpolation
// 	2	use bilinear interpolation
// 	3	use bicubic interpolation
//
// 	p	toggle periodic extrapolation
// 	w	toggle horizon
// 	.	toggle grid points
//
//
// Compilation:
//
//	cc -O3 -DNDEBUG viho.c iio.c ftr.c -lX11 -lpng -ljpeg -ltiff -o viho 


// SECTION 1. Libraries and data structures                                 {{{1

// standard libraries
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

//// user interface library
#include "ftr.h"


// radius of the disks that are displayed around control points
#define DISK_RADIUS 7

// zoom factor for zoom-in and zoom-out
#define ZOOM_FACTOR 1.43


// data structure to store the state of the viewer
struct viewer_state {
	// image data
	float *img;
	int iw, ih, pd;

	// geometry
	double c[4][2]; // control points in window coordinates
	double p[4][2]; // control points in image coordinates
	double offset[2], scale; // window viewport

	// dragging state
	bool dragging_window_point;
	bool dragging_image_point;
	bool dragging_background;
	int dragged_point;
	int drag_handle[2];

	// display options
	int interpolation_order; // 0=nearest, 1=linear, 2=bilinear, 3=bicubic
	bool tile_plane;
	bool show_horizon;
	bool show_grid_points;
	bool restrict_to_affine;
};


// function to reset and center the viewer
static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// set control points near the corners of the window
	double margin = 33;
	e->c[0][0] = margin;
	e->c[0][1] = margin;
	e->c[1][0] = f->w - margin;
	e->c[1][1] = margin;
	e->c[2][0] = f->w - margin;
	e->c[2][1] = f->h - margin;
	e->c[3][0] = margin;
	e->c[3][1] = f->h - margin;

	// set control points exactly at the corners of the image
	e->p[0][0] = 0;
	e->p[0][1] = 0;
	e->p[1][0] = e->iw - 1;
	e->p[1][1] = 0;
	e->p[2][0] = e->iw - 1;
	e->p[2][1] = e->ih - 1;
	e->p[3][0] = 0;
	e->p[3][1] = e->ih - 1;

	// drag state
	e->dragging_window_point = false;
	e->dragging_image_point = false;
	e->dragged_point = -1;
	e->dragging_background = false;

	// viewport
	e->offset[0] = 0;
	e->offset[1] = 0;
	e->scale = 1;

	// visualization options
	e->interpolation_order = 0;
	e->tile_plane = false;
	e->show_horizon = false;
	e->show_grid_points = false;
	e->restrict_to_affine = false;

	f->changed = 1;
}


// funtion to test whether a point is inside the window
static int insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}



// SECTION 2. Linear algebra                                                {{{1

// y = H(x)
static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

// compute the inverse homography (inverse of a 3x3 matrix)
static double invert_homography(double invH[3][3], double H[3][3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = ( a[4] * a[8] - a[5] * a[7] ) / det;
	r[1] = ( a[2] * a[7] - a[1] * a[8] ) / det;
	r[2] = ( a[1] * a[5] - a[2] * a[4] ) / det;
	r[3] = ( a[5] * a[6] - a[3] * a[8] ) / det;
	r[4] = ( a[0] * a[8] - a[2] * a[6] ) / det;
	r[5] = ( a[2] * a[3] - a[0] * a[5] ) / det;
	r[6] = ( a[3] * a[7] - a[4] * a[6] ) / det;
	r[7] = ( a[1] * a[6] - a[0] * a[7] ) / det;
	r[8] = ( a[0] * a[4] - a[1] * a[3] ) / det;
	return det;
}

// C = AoB, composition of two homographies (product of 3x3 matrices)
static void compose_homographies(double C[3][3], double A[3][3], double B[3][3])
{
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 3; k++)
			C[i][j] += A[i][k] * B[k][j];
	}
}

// Find the homography that changes the canonical projective basis
// into the given four points (x, y, z, w)
static void homography_from_four_points(double H[3][3],
		double x[2], double w[2], double z[2], double y[2])
{
	// fix the degree of freedom (assuming the four points are finite)
	double t = 1;

	// translation coefficients
	double p = x[0];
	double q = x[1];

	// "core" 2x2 system
	double A = w[0] - z[0];
	double B = y[0] - z[0];
	double C = w[1] - z[1];
	double D = y[1] - z[1];
	double P = z[0] - y[0] - w[0] + p;
	double Q = z[1] - y[1] - w[1] + q;
	double DET = A * D - B * C;
	double r = (D * P - B * Q) / DET;
	double s = (A * Q - C * P) / DET;

	// solve the rest of the diagonal system
	double a = w[0] * ( 1 + r ) - p;
	double b = y[0] * ( 1 + s ) - p;
	double c = w[1] * ( 1 + r ) - q;
	double d = y[1] * ( 1 + s ) - q;

	// fill-in the output
	H[0][0] = a; H[0][1] = b; H[0][2] = p;
	H[1][0] = c; H[1][1] = d; H[1][2] = q;
	H[2][0] = r; H[2][1] = s; H[2][2] = t;
}

// Find the homography that moves the four points (x,y,z,w) to (a,b,c,d)
static void homography_from_eight_points(double H[3][3],
		double x[2], double y[2], double z[2], double w[2],
		double a[2], double b[2], double c[2], double d[2])
{
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}

/// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}



// SECTION 3. Coordinate Conversions                                        {{{1

// This program deals with three systems of coordinates
//
// "image"  : coordinates in the rectangular domain of the input image
// "view"   : coordinates in the infinite plane where the image is mapped
// "window" : coordinates in the window, which is a rectangluar piece of "view"
//

// change from view coordinates to window coordinates
static void map_view_to_window(struct viewer_state *e, double y[2], double x[2])
{
	for (int k = 0; k < 2; k++)
		y[k] = e->offset[k] + e->scale * x[k];
}

// change from window coordinates to view coordinates
static void map_window_to_view(struct viewer_state *e, double y[2], double x[2])
{
	for (int k = 0; k < 2; k++)
		y[k] = ( x[k] - e->offset[k] ) / e->scale;
}

// obtain the direct homography from the current configuration
static void obtain_current_homography(double H[3][3], struct viewer_state *e)
{
	double C[4][2];
	for (int p = 0; p < 4; p++)
		map_view_to_window(e, C[p], e->c[p]);
	homography_from_eight_points(H,
			C[0], C[1], C[2], C[3],
			e->p[0], e->p[1], e->p[2], e->p[3]
			);
}

// change from window coordinates to image coordinates
static void map_window_to_image(struct viewer_state *e, double *y, double *x)
{
	double H[3][3], C[4][2];
	obtain_current_homography(H, e);
	apply_homography(y, H, x);
}

// Convert a floating-point color into a byte in the range [0,255]
// (Colors outisde of this range are saturated.)
static double float_to_byte(double x)
{
	if (x < 0) return 0;
	if (x < 255) return x;
	return 255;
}

static void adjust_affine_corner(struct viewer_state *e)
{
	for (int i = 0; i < 2; i++)
		e->c[2][i] = e->c[1][i] + ( e->c[3][i] - e->c[0][i] );
}


// SECTION 4. Extrapolation                                                 {{{1

// A "extrapolator" evaluates an image at an arbitrary integral position.
// When the position is outside the image domain, the value is extrapolated
// (by periodization, or by a constant value).

// type of the "extrapolator" functions
typedef float (*extrapolator_t)(float*,int,int,int,int,int,int);

// auxiliary function: compute n%p correctly, even for huge and negative numbers
static int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}

// instance of "extrapolator_t", extrapolate by periodicity
static float getsample_per(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// instance of "extrapolator_t", extrapolate by a constant value
static float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return value;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// obtain an extrapolator function from the current configuration
static extrapolator_t obtain_extrapolator(struct viewer_state *e)
{
	if (e->tile_plane)
		return getsample_per;
	float top = 255; // white
	getsample_cons(&top, 0, 0, 0, 0, 0, 0);
	return getsample_cons;
}



// SECTION 5. Local Interpolation                                           {{{1
//
// An "interpolator" evaluates an image at an arbitrary floating-point
// position.  When the position is not integral, the value is a combination of
// neighboring values.  Notice that when the position is outside the image
// domain, an extrapolator is also needed.

// type of the "interpolator" functions
typedef float (*interpolator_t)(float*,int,int,int,float,float,int,
		extrapolator_t);

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// instance of "interpolator_t", for bilinear interpolation
static float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}

// auxiliary code for "linear" interpolation
#include "marching_interpolation.c"

// instance of "interpolator_t", for linear (marching) interpolation
static float linear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return marchi(a, c, b, d, p-ip, q-iq);
}

// instance of "interpolator_t" for nearest neighbor interpolation
static float nearest_neighbor_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = round(p);
	int iq = round(q);
	return pix(x, w, h, pd, ip, iq, l);
}


// one-dimensional cubic interpolation of four data points ("Keys")
static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

// two-dimensional separable cubic interpolation, on a 4x4 grid
static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

// instance of "interpolator_t" for bicubic interpolation
static float bicubic_interpolation_at(float *img, int w, int h, int pd,
		float x, float y, int l, extrapolator_t p)
{
	x -= 1;
	y -= 1;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, pd, ix + i, iy + j, l);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}

// obtain an interpolation operator from the current configuration
static interpolator_t obtain_interpolator(struct viewer_state *e)
{
	if (e->dragged_point >= 0 || e->dragging_background)
		return nearest_neighbor_at;
	if (e->interpolation_order == 1) return linear_interpolation_at;
	if (e->interpolation_order == 2) return bilinear_interpolation_at;
	if (e->interpolation_order == 3) return bicubic_interpolation_at;
	return nearest_neighbor_at;
}



// SECTION 6. Main Warping Function                                         {{{1

// draw the image warped by the current homography
static void draw_warped_image(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	int w = e->iw;
	int h = e->ih;

	double         H[3][3];   obtain_current_homography(H, e);
	extrapolator_t OUT      = obtain_extrapolator(e);
	interpolator_t EVAL     = obtain_interpolator(e);

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2] = {i, j};
		apply_homography(p, H, p);
		p[0] = p[0] * w / (w - 1.0) - 0.5;
		p[1] = p[1] * h / (h - 1.0) - 0.5;
		for (int l = 0; l < 3; l++)
		{
			int idx = l + 3 * (f->w * j + i);
			float v = EVAL(e->img, w, h, e->pd, p[0], p[1], l, OUT);
			f->rgb[idx] = float_to_byte(v);
		}
	}
}



// SECTION 7. Drawing                                                       {{{1

// Subsection 7.1. Drawing segments                                         {{{2

// generic function to traverse a segment between two pixels
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (qx - px > qy - py || px - qx > qy - py) { // horizontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px + j*slope), j+py, e);
		}
	}
}

// auxiliary function for drawing a red pixel
static void plot_pixel_red(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 255;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
	}
}

// auxiliary function for drawing a green pixel
static void plot_pixel_green(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 128;
		f->rgb[3*idx+2] = 0;
	}
}

// function to draw a red segment
static void plot_segment_red(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

// function to draw a green segment
static void plot_segment_green(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_green, f);
}



// Subsection 7.2. Drawing user-interface elements                          {{{2

// draw the positions of the image samples
static void draw_grid_points(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	double H[3][3], iH[3][3];
	obtain_current_homography(iH, e);
	invert_homography(H, iH);
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 0; // black
	for (int j = 0; j < e->ih; j++)
	for (int i = 0; i < e->iw; i++)
	{
		double p[2] = {i, j};
		apply_homography(p, H, p);
		int ip[2] = {p[0], p[1]};
		if (insideP(f, ip[0], ip[1]))
		for (int l = 0; l < 3; l++)
		{
			int idx_win = l + 3 * (f->w * ip[1] + ip[0]);
			float val = getsample_cons(e->img, e->iw, e->ih,
					e->pd, i, j, l);
			f->rgb[idx_win] = 127+val/2;
		}
	}
}

// draw four red segments connecting the control points
static void draw_four_red_segments(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	for (int p = 0; p < 4; p++)
	{
		double x0[2], xf[2];
		map_view_to_window(e, x0, e->c[p]);
		map_view_to_window(e, xf, e->c[(p+1)%4]);
		plot_segment_red(f, x0[0], x0[1], xf[0], xf[1]);
	}
}

// draw disks around the control points
static void draw_four_control_points(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	for (int p = 0; p < 4; p++)
	{
		if (p == 2 && e->restrict_to_affine) continue;
		double P[2];
		map_view_to_window(e, P, e->c[p]);

		// grey circle
		int side = DISK_RADIUS;
		for (int j = -side-1 ; j <= side+1; j++)
		for (int i = -side-1 ; i <= side+1; i++)
		if (hypot(i, j) < side)
		{
			int ii = P[0] + i;
			int jj = P[1] + j;
			if (insideP(f, ii, jj))
				for (int c = 0; c < 3; c++)
					f->rgb[3*(f->w*jj+ii)+c] = 127;
		}

		// central green dot
		int ii = P[0];
		int jj = P[1];
		if (insideP(f, ii, jj))
			f->rgb[3*(f->w*jj+ii)+1]=255;
	}
}

// plot the horizon in green
static void draw_horizon(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// four control points in projective coordinates
	double a[3] = {e->c[0][0], e->c[0][1], 1};
	double b[3] = {e->c[1][0], e->c[1][1], 1};
	double c[3] = {e->c[2][0], e->c[2][1], 1};
	double d[3] = {e->c[3][0], e->c[3][1], 1};

	// lines though the control points
	double ab[3]; vector_product(ab, a, b);
	double cd[3]; vector_product(cd, c, d);
	double ad[3]; vector_product(ad, a, d);
	double bc[3]; vector_product(bc, b, c);

	// intersections of opposite sides (vanishing points)
	double p[3]; vector_product(p, ab, cd);
	double q[3]; vector_product(q, ad, bc);

	// horizon := line through two vanishing points
	double horizon[3]; vector_product(horizon, p, q);

	// affine coordinates of points
	for (int k = 0; k < 2; k++) {
		p[k] /= p[2];
		q[k] /= q[2];
	}

	// plot the horizon
	double v[2][2];
	map_view_to_window(e, v[0], p);
	map_view_to_window(e, v[1], q);
	if (hypot(hypot(v[0][0], v[0][1]), hypot(v[1][0], v[1][1])) < 1e5)
		plot_segment_green(f, v[0][0], v[0][1], v[1][0], v[1][1]);

}

// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	if (e->show_grid_points)
		draw_grid_points(f);
	else
		draw_warped_image(f);
	draw_four_red_segments(f);
	draw_four_control_points(f);
	if (e->show_horizon)
		draw_horizon(f);
}



// SECTION 8. User-Interface Actions and Events                             {{{1

// action: viewport translation
static void change_view_offset(struct viewer_state *e, double dx, double dy)
{
	e->offset[0] += dx;
	e->offset[1] += dy;
}

// action: viewport zoom
static void change_view_scale(struct viewer_state *e, int x, int y, double fac)
{
	double center[2], X[2] = {x, y};
	map_window_to_view(e, center, X);
	e->scale *= fac;
	for (int p = 0; p < 2; p++)
		e->offset[p] = -center[p]*e->scale + X[p];
	fprintf(stderr, "zoom changed %g\n", e->scale);
}

// action: drag a point in the window domain
static void drag_point_in_window_domain(struct viewer_state *e, int x, int y)
{
	int p = e->dragged_point;
	double X[2] = {x, y};
	map_window_to_view(e, e->c[p], X);
	if (e->restrict_to_affine)
		adjust_affine_corner(e);
}

// action: drag a point in the image domain
static void drag_point_in_image_domain(struct viewer_state *e, int x, int y)
{
	int p = e->dragged_point;
	double P[2], Q[2] = {x, y};
	map_window_to_image(e, P, Q);
	e->p[p][0] = P[0];
	e->p[p][1] = P[1];
	map_window_to_view(e, e->c[p], Q);
}

// test whether (x,y) is inside one of the four control disks
static int hit_point(struct viewer_state *e, double x, double y)
{
	for (int p = 0; p < 4; p++)
	{
		if (p == 2 && e->restrict_to_affine) continue;
		double P[2];
		map_view_to_window(e, P, e->c[p]);
		if (hypot(P[0] - x, P[1] - y) < 2 + DISK_RADIUS)
			return p;
	}
	return -1;
}

// key handler
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q') {
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'c') center_view(f);
	if (k == 'j') change_view_offset(e, 0, -10);
	if (k == 'k') change_view_offset(e, 0, 10);
	if (k == 'h') change_view_offset(e, 10, 0);
	if (k == 'l') change_view_offset(e, -10, 0);
	if (k == FTR_KEY_DOWN ) change_view_offset(e, 0, -100);
	if (k == FTR_KEY_UP   ) change_view_offset(e, 0, 100);
	if (k == FTR_KEY_RIGHT) change_view_offset(e, -100, 0);
	if (k == FTR_KEY_LEFT)  change_view_offset(e, 100, 0);
	if (k == '+') change_view_scale(e, f->w/2, f->h/2, ZOOM_FACTOR);
	if (k == '-') change_view_scale(e, f->w/2, f->h/2, 1.0/ZOOM_FACTOR);
	if (k == 'p') e->tile_plane = !e->tile_plane;
	if (k == 'w') e->show_horizon = !e->show_horizon;
	if (k >= '0' && k <= '9') e->interpolation_order = k - '0';
	if (k == '.') e->show_grid_points = !e->show_grid_points;
	if (k == 'a') e->restrict_to_affine = !e->restrict_to_affine;

	e->dragging_window_point = e->dragging_image_point = false;
	f->changed = 1;
}

// resize handler
static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	f->changed = 1;
}

// mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	int p = hit_point(e, x, y);

	// begin dragging a control point in the WINDOW DOMAIN
	if (k == FTR_BUTTON_LEFT && p >= 0)
	{
		e->dragged_point = p;
		e->dragging_window_point = true;
	}

	// end dragging a control point in the WINDOW DOMAIN
	if (e->dragging_window_point && k == -FTR_BUTTON_LEFT)
	{
		drag_point_in_window_domain(e, x, y);
		e->dragging_window_point = false;
		e->dragged_point = -1;
	}

	// begin dragging a control point in the IMAGE DOMAIN
	if (k == FTR_BUTTON_RIGHT && p >= 0)
	{
		e->dragged_point = p;
		e->dragging_image_point = true;
	}

	// end dragging a control point in the IMAGE DOMAIN
	if (e->dragging_image_point && k == -FTR_BUTTON_RIGHT)
	{
		drag_point_in_image_domain(e, x, y);
		e->dragging_image_point = false;
		e->dragged_point = -1;
	}

	// begin dragging a the WINDOW BACKGROUND
	if (k == FTR_BUTTON_LEFT && hit_point(e, x, y) < 0)
	{
		e->drag_handle[0] = x;
		e->drag_handle[1] = y;
		e->dragging_background = true;
	}

	// end dragging the WINDOW BACLGROUND
	if (e->dragging_background && k == -FTR_BUTTON_LEFT)
	{
		int dx = x - e->drag_handle[0];
		int dy = y - e->drag_handle[1];
		change_view_offset(e, dx, dy);
		e->dragging_background = false;
	}

	// zoom in/out
	if (k == FTR_BUTTON_DOWN) change_view_scale(e, x, y, ZOOM_FACTOR);
	if (k == FTR_BUTTON_UP) change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);

	f->changed = 1;
}

// mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// drag WINDOW DOMAIN control point (realtime feedback)
	if (e->dragging_window_point && m & FTR_BUTTON_LEFT)
	{
		drag_point_in_window_domain(e, x, y);
		f->changed = 1;
	}

	// drag IMAGE DOMAIN control point (realtime feedback)
	if (e->dragging_image_point && m & FTR_BUTTON_RIGHT)
	{
		drag_point_in_image_domain(e, x, y);
		f->changed = 1;
	}

	// drag WINDOW DOMAIN background (realtime feedback)
	if (e->dragging_background && m & FTR_BUTTON_LEFT)
	{
		int dx = x - e->drag_handle[0];
		int dy = y - e->drag_handle[1];
		change_view_offset(e, dx, dy);
		e->drag_handle[0] = x;
		e->drag_handle[1] = y;
		f->changed = 1;
	}
}

// expose handler
static void event_expose(struct FTR *f, int b, int m, int x, int y)
{
	if (f->changed)
		paint_state(f);
}


// SECTION 9. Image processing

//#include "xmalloc.c"
//
//#define PYRAMID_LEVELS 20
//
//struct image_pyramid {
//	float *x[PYRAMID_LEVELS];
//	int w[PYRAMID_LEVELS];
//	int h[PYRAMID_LEVELS];
//	int pd;
//};
//
//static void do_pyramid(struct image_pyramid *p, float *x, int w, int h, int pd)
//{
//	fprintf(stderr, "building a multiscale pyramid %d x %d x %d\n",
//			w, h, pd);
//
//	int nw = 1 + ceil(log2(w+1));
//	int nh = 1 + ceil(log2(h+1));
//	p->n = (nw + nh + abs(nw - nh))/2;
//	p->x = xmalloc(n * sizeof*p->x);
//	p->w = xmalloc(n * sizeof*p->w);
//	p->h = xmalloc(n * sizeof*p->h);
//
//	p->w[0] = w;
//	p->h[0] = h;
//	for (int i = 1; i < p->n; i++)
//
//
//	fprintf(stderr, "\t%d x %d\n", nw, nh);
//}
//
//static void free_pyramid(struct image_pyramid *p)
//{
//	for (int i = 0; i < pyra i++)
//		free(p->x[i]);
//	free(p->x);
//	free(
//}



// SECTION 10. Main Program                                                 {{{1

// library for image input-output
#include "iio.h"

// main function
int main(int argc, char *argv[])
{
	if (argc != 2 && argc != 1) {
		fprintf(stderr, "usage:\n\t%s [image.png]\n", *argv);
		return 1;
	}
	char *filename_in = argc == 2 ? argv[1] : "/tmp/lenak.png";

	struct FTR f = ftr_new_window(512,512);
	struct viewer_state e[1];
	f.userdata = e;

	e->img = iio_read_image_float_vec(filename_in, &e->iw, &e->ih, &e->pd);

	//struct image_pyramid p[1];
	//do_pyramid(p, e->img, e->iw, e->ih, e->pd);

	center_view(&f);

	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "resize", event_resize);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "key", event_key);

	return ftr_loop_run(&f);
}

// vim:set foldmethod=marker:
