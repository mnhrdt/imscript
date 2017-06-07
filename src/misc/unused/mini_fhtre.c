#include <stdio.h>
#include <stdlib.h>
//#include "library_for_reading_and_writing_images.h"
#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 4)
		return printf("usage:\n\t%s threshold in out\n", *argv);

	int width, height;
	float threshold = atof(argv[1]);
	float *image = read_image_into_floats(argv[2], &width, &height);
	for (int i = 0; i < width*height; i++)
		image[i] = image[i] > threshold;
	write_floats_into_image(argv[3], image, width, height);

	return 0;
}
