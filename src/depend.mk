plambda.o: plambda.c smapa.h fail.c xmalloc.c random.c parsenumbers.c \
 colorcoordsf.c getpixel.c iio.h
vecov.o: vecov.c iio.h fail.c xmalloc.c random.c smapa.h help_stuff.c \
 pickopt.c
veco.o: veco.c iio.h fail.c xmalloc.c random.c help_stuff.c pickopt.c
vecoh.o: vecoh.c iio.h fail.c xmalloc.c random.c modes_detector.c smapa.h \
 help_stuff.c pickopt.c
morsi.o: morsi.c iio.h help_stuff.c
downsa.o: downsa.c iio.h fail.c xmalloc.c random.c help_stuff.c
upsa.o: upsa.c iio.h fail.c marching_squares.c marching_interpolation.c \
 bicubic.c getpixel.c
ntiply.o: ntiply.c iio.h
censust.o: censust.c iio.h pickopt.c
dither.o: dither.c iio.h
qauto.o: qauto.c iio.h help_stuff.c pickopt.c
qeasy.o: qeasy.c iio.h
homwarp.o: homwarp.c extrapolators.c bilinear_interpolation.c \
 marching_interpolation.c bicubic_gray.c spline.c iio.h xmalloc.c fail.c \
 parsenumbers.c help_stuff.c pickopt.c
synflow.o: synflow.c iio.h xmalloc.c fail.c synflow_core.c getpixel.c \
 marching_interpolation.c vvector.h homographies.c smapa.h
backflow.o: backflow.c iio.h fail.c xmalloc.c getpixel.c bicubic.c \
 smapa.h
flowinv.o: flowinv.c iio.h fail.c xmalloc.c bicubic.c getpixel.c
nnint.o: nnint.c abstract_heap.h xmalloc.c fail.c iio.h pickopt.c
bdint.o: bdint.c abstract_dsf.c iio.h pickopt.c
amle.o: amle.c iio.h fail.c xmalloc.c smapa.h
simpois.o: simpois.c cleant_cgpois.c minicg.c smapa.h iio.h pickopt.c
ghisto.o: ghisto.c iio.h xmalloc.c fail.c smapa.h
contihist.o: contihist.c xfopen.c fail.c xmalloc.c iio.h
fontu.o: fontu.c xmalloc.c fail.c xfopen.c dataconv.c iio.h pickopt.c
imprintf.o: imprintf.c iio.h help_stuff.c
pview.o: pview.c iio.h fail.c xmalloc.c xfopen.c parsenumbers.c \
 drawsegment.c pickopt.c smapa.h random.c
viewflow.o: viewflow.c iio.h fail.c drawsegment.c colorcoordsf.c \
 marching_squares.c
flowarrows.o: flowarrows.c iio.h fail.c xmalloc.c drawsegment.c \
 getpixel.c smapa.h
palette.o: palette.c fail.c xmalloc.c xfopen.c smapa.h iio.h
ransac.o: ransac.c fail.c xmalloc.c xfopen.c random.c ransac_cases.c \
 vvector.h homographies.c moistiv_epipolar.c parsenumbers.c
srmatch.o: srmatch.c fail.c xmalloc.c xfopen.c siftie.c parsenumbers.c \
 smapa.h ok_list.c grid.c iio.h ransac.c random.c ransac_cases.c \
 vvector.h homographies.c moistiv_epipolar.c
tiffu.o: tiffu.c
siftu.o: siftu.c siftie.c fail.c xmalloc.c xfopen.c parsenumbers.c \
 smapa.h ok_list.c grid.c iio.h
crop.o: crop.c fail.c xmalloc.c iio.h
lrcat.o: lrcat.c iio.h xmalloc.c fail.c getpixel.c pickopt.c
tbcat.o: tbcat.c iio.h smapa.h xmalloc.c fail.c getpixel.c pickopt.c
fftshift.o: fftshift.c iio.h
imflip.o: imflip.c help_stuff.c iio.h
bmms.o: bmms.c xmalloc.c fail.c getpixel.c iio.h pickopt.c
registration.o: registration.c iio.h
blur.o: blur.c fail.c xmalloc.c smapa.h iio.h parsenumbers.c
fft.o: fft.c iio.h fail.c xmalloc.c
dct.o: dct.c iio.h
dht.o: dht.c iio.h xmalloc.c fail.c
flambda.o: flambda.c smapa.h fail.c xmalloc.c random.c parsenumbers.c \
 colorcoordsf.c fancy_image.h getpixel.c
fancy_crop.o: fancy_crop.c fancy_image.h
fancy_downsa.o: fancy_downsa.c fancy_image.h
iion.o: iion.c iio.h
