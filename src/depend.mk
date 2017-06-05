plambda: plambda.c smapa.h fail.c xmalloc.c random.c parsenumbers.c \
 colorcoordsf.c getpixel.c iio.h
vecov: vecov.c iio.h fail.c xmalloc.c random.c smapa.h help_stuff.c \
 pickopt.c
veco: veco.c iio.h fail.c xmalloc.c random.c help_stuff.c pickopt.c
vecoh: vecoh.c iio.h fail.c xmalloc.c random.c modes_detector.c smapa.h \
 help_stuff.c pickopt.c
morsi: morsi.c iio.h help_stuff.c
downsa: downsa.c iio.h fail.c xmalloc.c random.c help_stuff.c
upsa: upsa.c iio.h fail.c marching_squares.c marching_interpolation.c \
 bicubic.c getpixel.c
ntiply: ntiply.c iio.h
censust: censust.c iio.h pickopt.c
dither: dither.c iio.h
qauto: qauto.c iio.h help_stuff.c pickopt.c
qeasy: qeasy.c iio.h
homwarp: homwarp.c extrapolators.c bilinear_interpolation.c \
 marching_interpolation.c bicubic_gray.c spline.c iio.h xmalloc.c fail.c \
 parsenumbers.c help_stuff.c pickopt.c
synflow: synflow.c iio.h xmalloc.c fail.c synflow_core.c getpixel.c \
 marching_interpolation.c vvector.h homographies.c smapa.h
backflow: backflow.c iio.h fail.c xmalloc.c getpixel.c bicubic.c \
 smapa.h
flowinv: flowinv.c iio.h fail.c xmalloc.c bicubic.c getpixel.c
nnint: nnint.c abstract_heap.h xmalloc.c fail.c iio.h pickopt.c
bdint: bdint.c abstract_dsf.c iio.h pickopt.c
amle: amle.c iio.h fail.c xmalloc.c smapa.h
simpois: simpois.c cleant_cgpois.c minicg.c smapa.h iio.h pickopt.c
ghisto: ghisto.c iio.h xmalloc.c fail.c smapa.h
contihist: contihist.c xfopen.c fail.c xmalloc.c iio.h
fontu: fontu.c xmalloc.c fail.c xfopen.c dataconv.c iio.h pickopt.c
imprintf: imprintf.c iio.h help_stuff.c
pview: pview.c iio.h fail.c xmalloc.c xfopen.c parsenumbers.c \
 drawsegment.c pickopt.c smapa.h random.c
viewflow: viewflow.c iio.h fail.c drawsegment.c colorcoordsf.c \
 marching_squares.c
flowarrows: flowarrows.c iio.h fail.c xmalloc.c drawsegment.c \
 getpixel.c smapa.h
palette: palette.c fail.c xmalloc.c xfopen.c smapa.h iio.h
ransac: ransac.c fail.c xmalloc.c xfopen.c random.c ransac_cases.c \
 vvector.h homographies.c moistiv_epipolar.c parsenumbers.c
srmatch: srmatch.c fail.c xmalloc.c xfopen.c siftie.c parsenumbers.c \
 smapa.h ok_list.c grid.c iio.h ransac.c random.c ransac_cases.c \
 vvector.h homographies.c moistiv_epipolar.c
tiffu: tiffu.c
siftu: siftu.c siftie.c fail.c xmalloc.c xfopen.c parsenumbers.c \
 smapa.h ok_list.c grid.c iio.h
crop: crop.c fail.c xmalloc.c iio.h
lrcat: lrcat.c iio.h xmalloc.c fail.c getpixel.c pickopt.c
tbcat: tbcat.c iio.h smapa.h xmalloc.c fail.c getpixel.c pickopt.c
fftshift: fftshift.c iio.h
imflip: imflip.c help_stuff.c iio.h
bmms: bmms.c xmalloc.c fail.c getpixel.c iio.h pickopt.c
registration: registration.c iio.h
blur: blur.c fail.c xmalloc.c smapa.h help_stuff.c parsenumbers.c iio.h
fft: fft.c iio.h fail.c xmalloc.c
dct: dct.c iio.h
dht: dht.c iio.h xmalloc.c fail.c
flambda: flambda.c smapa.h fail.c xmalloc.c random.c parsenumbers.c \
 colorcoordsf.c fancy_image.h getpixel.c
fancy_crop: fancy_crop.c fancy_image.h
fancy_downsa: fancy_downsa.c fancy_image.h
iion: iion.c iio.h
iion_u16: iion_u16.c iio.h
