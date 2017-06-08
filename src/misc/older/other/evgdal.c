// gcc geotif_getpixel.c iio.o -lgdal -ltiff -lpng -ljpeg
#include <stdio.h>
#include <gdal/gdal.h>
#include "fancy_image.c"

int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr,
			"usage:\n\t%s geotiff_image.tif < points.txt\n", *v);

	// get geotif tags with gdal (it's ugly but geotif API documentation sucks)
	GDALAllRegister();
	GDALDatasetH gdal_dataset = GDALOpen(v[1], GA_ReadOnly);
	if (gdal_dataset == NULL) {
		fprintf(stderr, "GDALOpen failed\n");
		return EXIT_FAILURE;
	}
	double tmp[6];
	double origin[2] = {0, 0};
	double scale[2] = {1, 1};
	if (GDALGetGeoTransform(gdal_dataset, tmp) == CE_None) {
		origin[0] = tmp[0], origin[1] = tmp[3];
		scale[0] = tmp[1], scale[1] = tmp[5];
	} else {
		fprintf(stderr, "cannot read origin and scale information\n");
	}

	// convert utm coordinates to pixels and read image samples
	struct fancy_image *a = fancy_image_open(v[1], "r");
	float utm[2], pix[2];
	while (2 == scanf("%g %g\n", utm, utm+1)) {
		pix[0] = (utm[0] - origin[0]) / scale[0];
		pix[1] = (utm[1] - origin[1]) / scale[1];
		//printf("pixel coordinates: %g %g\n", pix[0], pix[1]);
		for (int i = 0; i < a->pd; i++)
			printf("%g%c",
				fancy_image_getsample(a, pix[0], pix[1], i),
					i == a->pd - 1 ? '\n' : ' ');
	}
	return EXIT_SUCCESS;
}
