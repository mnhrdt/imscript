#include <math.h>

// returns a distance, in meters, to the center of mass of the earth
double wgs84(double lat, double lon)
{
	double a = 6378137;        // semi-major axis
	double f = 298.257223563;  // inverse-flatting
	double b = a*(1- 1/f);     // semi-minor axis
}


#ifdef MAIN_WGS84
#include <stdio.h>
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s longitude latitude\n", *v);
		return 1;
	}
	double lon = atof(v[1]);
	double lat = atof(v[2]);
	double r = wgs84(lon, lat);
	printf("%g\n", r);
	return 0;
}
#endif//MAIN_SRTM4
