#include <math.h>

#include "xmalloc.c"


// out[0] = K
// out[1] = G^x_xx
// out[2] = G^x_xy
// out[3] = G^x_yy
// out[4] = G^y_xx
// out[5] = G^y_xy
// out[6] = G^y_yy
void darboux(float *kchrist, float *EFG, int w, int h)
{
	//int n = w*h*sizeof(float);
	//float *Ex = xmalloc(n); float *Ey = xmalloc(n);
	//float *Fx = xmalloc(n); float *Fy = xmalloc(n);
	//float *Gx = xmalloc(n); float *Gy = xmalloc(n);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float E = p(EFG,w,h,3, i,j,0);
		float F = p(EFG,w,h,3, i,j,1);
		float G = p(EFG,w,h,3, i,j,2);
		float Ex = p(EFG,w,h,3, i+1,j,0)-E;
		float Ey = p(EFG,w,h,3, i,j+1,0)-E;
		float Fx = p(EFG,w,h,3, i+1,j,1)-F;
		float Fy = p(EFG,w,h,3, i,j+1,1)-F;
		float Gx = p(EFG,w,h,3, i+1,j,2)-G;
		float Gy = p(EFG,w,h,3, i,j+1,2)-G;
		float disc = 2*(E*G-F*F);
		float G_xxx = ( G*Ex - 2*F*Fx + F*Ey )/disc;
		float G_xxy = ( G*Ey - F*Gx          )/disc;
		float G_xyy = ( 2*G*Fy - G*Gx - F*Gy )/disc;
		float G_yxx = ( 2*E*Fx - E*Ey - F*Ex )/disc;
		float G_yxy = ( E*Gx - F*Ey          )/disc;
		float G_xxx = ( E*Gy - 2*F*Fy + F*Gx )/disc;
		float Evv = p(EFG,w,h,3, i,j-1,0)-2*E+p(EFG,w,h,3, i,j+1,0);
		float Gvv = p(EFG,w,h,3, i,j-1,0)-2*E+p(EFG,w,h,3, i,j+1,0);
		float K = NAN;
	}
	//free(Ex); free(Ey); free(Fx); free(Fy); free(Gx); free(Gy);
}

#include <stdio.h>

#include "iio.h"
int main(int c, char *v)
{
	if (c != 1 && c != 2 && c != 3)
	{
		fprintf(stderr, "usage:\n\t%s [metric [christoffel]]\n", *v);
		//                          0  1       2
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *EFG = xmalloc(w*h*3*sizeof*EFG);
	get_EFG(EFG, x, w, h, pd);
	float *kchrist = xmalloc(w*h*7*sizeof*kchrist);
}
