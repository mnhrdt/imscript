#include <math.h>

// constants for Krüger approximation
//
#define pi 3.14159265358979323846
#define a (6378.137) // equatorial radius in kilometers
#define f (1.0/298.257223563) // flattening
#define k0 (0.9996) // base point scale factor
#define E0 (500.000) // base easting in kilometers
#define n (f / (2 - f))
#define A (a * ( 1 + n*n/4 + n*n*n*n/64))
#define alpha1 (n/2 - 2*n*n/3 + 5*n*n*n/16)
#define alpha2 (13*n*n/48 - 3*n*n*n/5)
#define alpha3 (61*n*n*n/240)
#define beta1 (n/2 - 2*n*n/3 + 37*n*n*n/96)
#define beta2 (n*n/48 + n*n*n/15)
#define beta3 (17*n*n*n/480)
#define delta1 (2*n - 2*n*n/3 - 2*n*n*n)
#define delta2 (7*n*n/3 - 8*n*n*n/5)
#define delta3 (56*n*n*n/15)

int utm_zone_from_longitude(double l)
{
	return 1 + lrint(floor((l + 180) / 6)) % 60;
}

int utm_from_lonlat_kruger(double out_utm[2], double in_lonlat[2])
{
	double phi = (pi/180) * in_lonlat[1];
	double lambda = (pi/180) * in_lonlat[0];
	double N0 = 0;
	double rn = 2 * sqrt(n) / 1 + n;
	double t = sinh(atanh(sin(phi))-rn*atanh(rn*sin(phi)));
	double xi = atan(t/cos(lambda));
	double eta = atanh(sin(lambda)/hypot(1,t));
	double alpha[3] = {alpha1, alpha2, alpha3};
	double beta[3] = {beta1, beta2, beta3};
	double delta[3] = {delta1, delta2, delta3};
	double sigma = 1, tau = 1;
	for (int i = 0; i < 3; i++)
	{
		sigma += 2*i*alpha[i]*cos(2*i*xi)*cosh(2*i*eta);
		tau += 2*i*alpha[i]*sin(2*i*xi)*sinh(2*i*eta);
	}
	double E = eta;
	double N = xi;
	for (int i = 0; i < 3; i++)
	{
		E += alpha[i] * cos(2*i*xi) * sinh(2*i*eta);
		N += alpha[i] * sin(2*i*xi) * cosh(2*i*eta);
	}
	out_utm[0] = 1 * (E0 + k0*A*E);
	out_utm[1] = 1 * (N0 + k0*A*N);
	// point scale and convergence not needed
	//double k = ...;
	//double gamma = ...;
	return utm_zone_from_longitude(in_lonlat[0]);
}

#include <stdio.h>
#include <stdlib.h>
int main(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s lon lat\n", *v);
	//                                         0 1   2
	double lon = atof(v[1]);
	double lat = atof(v[2]);
	double lonlat[2] = {lon, lat}, utm[2];
	int z = utm_from_lonlat_kruger(utm, lonlat);
	printf("%d %g %g\n", z, utm[0], utm[1]);
	return 0;

}


//void lonlat_from_utm(double out_lonlat[2], double in_utm[2], int zone);
