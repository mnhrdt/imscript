// gcc geotiff_getpixel.c iio.o -lgeotiff -ltiff -lpng -ljpeg
#include <stdio.h>
#include <geo_tiffp.h>
#include "fancy_image.c"

int main(int c, char *v[])
{
    if (c != 2)
        return fprintf(stderr, "usage:\n\t%s geotiff_image.tif < points.txt\n", *v);

    // get TIFF tags
    TIFF *tif = XTIFFOpen(v[1], "r");
    if (!tif)
        fail("failed in XTIFFOpen\n");

    int count;
    double *tags = NULL;
    double scale[2] = {1, -1};
    if (TIFFGetField(tif, GTIFF_PIXELSCALE, &count, &tags)) {
        scale[0] = tags[0];
        scale[1] = tags[1];
        fprintf(stderr, "scale %g, %g\n", scale[0], scale[1]);
    } else
        fprintf(stderr, "GTIFF_PIXELSCALE not defined\n");

    double origin[2] = {0, 0};
    if (TIFFGetField(tif, GTIFF_TIEPOINTS, &count, &tags)) {
        origin[0] = tags[3];
        origin[1] = tags[4];
        fprintf(stderr, "origin %g, %g\n", origin[0], origin[1]);
    } else
        fprintf(stderr, "GTIFF_TIEPOINTS not defined\n");

    XTIFFClose(tif);

    // convert utm coordinates to pixels and read image samples
    struct fancy_image *a = fancy_image_open(v[1], "r,megabytes=800");
    float utm[2], pix[2];
    while (2 == scanf("%g %g\n", utm, utm+1)) {
        pix[0] = (utm[0] - origin[0]) / scale[0];
        pix[1] = -(utm[1] - origin[1]) / scale[1];
        //printf("pixel coordinates: %g %g\n", pix[0], pix[1]);
        for (int i = 0; i < a->pd; i++)
            printf("%g%c", fancy_image_getsample(a, pix[0], pix[1], i),
                   i == a->pd - 1 ? '\n' : ' ');
    }
    return EXIT_SUCCESS;
}
