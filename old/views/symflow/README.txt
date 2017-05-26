Some standalone programs for manipulating optical flows.

********
BUILDING
********

Requirements: ANSI C compiler, libtiff, libjpeg, libpng

The scripts "build.sh", "test.sh" and "clean.sh" allow you to build, test and
clean the executables.


*****
TOOLS
*****


To see a mini-documentation for each tool, run the program without any
argument.

rgfieldst: generate a sequence of correlated random vector fields
viewflow: display a vector field in colors
flowarrows: display a vector field using arrows
backflow: warp an image by a vector field
flowinv: find the inverse of a vector field
qeasy: quantize an image within the specified dynamic range
plambda: evaluate an algebraic expression on images



*******
EXAMPLE
*******

The following unix commands simulate a turbulent vide of 100 frames with base
image lena.png (whose size is 256x256).

./rgfieldst 256 256 100  50    40 0 0    80 -70   80 field_%03d.flo
for i in {000..099}; do
	./backflow field_$i.flo lena.png | qeasy 0 255 - turbulent_lena_$i.png
done

The meaning of the arguments in the call to "rgfieldst" is the following:

256 256 = dimensions in pixels of each output image
100 = number of output images
50 = strength of the turbulence (unnormalized)
40 0 0, 80 -70, 80 = covariance matrix of the turbulence
field_%03d.flo = printf-like pattern for the filenames of the output images

With these parameters we obtain an upwards-moving pattern of turbulence with a
circular grain of size 10 and max strength about 3 pixels.  The shape and size
of the grain, and the direction and speed of the turbulent pattern are defined
by the coefficients of the covariance matrix.  This matrix must be symmetric
and positive definite.  The 3D ellipsoid determined by this matrix contains the
necessary information.  The ellipse of intersection of this ellipsoid with the
(x,y) plane determines the size and shape of the grain.  The spatiotemporal
direction of the rest of the ellipsoid determines the direction and speed of
movement of the pattern.  The length of the ellipsoid in the spatiotemporal
direction determines the speed at which the pattern may change in time.

Setting the coefficients of the covariance matrix is general enough but
admittedly non-intuitive.  By user request, I can provide a different and
easier interface.

Enric Meinardt-Llopis
