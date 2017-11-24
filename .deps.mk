src/abstract_dsf.o: src/abstract_dsf.c
src/amle.o: src/amle.c src/iio.h src/fail.c src/xmalloc.c src/smapa.h
src/autotrim.o: src/autotrim.c src/iio.h
src/backflow.o: src/backflow.c src/iio.h src/fail.c src/xmalloc.c \
 src/getpixel.c src/bicubic.c src/smapa.h
src/bdint.o: src/bdint.c src/abstract_dsf.c src/iio.h src/pickopt.c
src/bicubic.o: src/bicubic.c src/getpixel.c
src/bicubic_gray.o: src/bicubic_gray.c
src/bilinear_interpolation.o: src/bilinear_interpolation.c
src/blur.o: src/blur.c src/fail.c src/xmalloc.c src/smapa.h src/help_stuff.c \
 src/parsenumbers.c src/pickopt.c src/iio.h
src/bmms.o: src/bmms.c src/xmalloc.c src/fail.c src/getpixel.c src/iio.h \
 src/pickopt.c
src/ccproc.o: src/ccproc.c src/abstract_dsf.c src/xmalloc.c src/fail.c
src/censust.o: src/censust.c src/iio.h src/pickopt.c
src/cleant_cgpois.o: src/cleant_cgpois.c src/minicg.c src/smapa.h
src/colorcoordsf.o: src/colorcoordsf.c
src/colormatch.o: src/colormatch.c src/iio.h
src/contihist.o: src/contihist.c src/xfopen.c src/fail.c src/xmalloc.c \
 src/iio.h src/pickopt.c
src/crop.o: src/crop.c src/fail.c src/xmalloc.c src/iio.h
src/dataconv.o: src/dataconv.c src/xmalloc.c src/fail.c
src/dct.o: src/dct.c src/iio.h
src/dht.o: src/dht.c src/iio.h src/xmalloc.c src/fail.c
src/dither.o: src/dither.c src/iio.h
src/downsa.o: src/downsa.c src/iio.h src/fail.c src/xmalloc.c src/random.c \
 src/help_stuff.c
src/drawsegment.o: src/drawsegment.c
src/drawtriangle.o: src/drawtriangle.c
src/exterior_algebra.o: src/exterior_algebra.c
src/extrapolators.o: src/extrapolators.c
src/fail.o: src/fail.c
src/fancy_crop.o: src/fancy_crop.c src/fancy_image.h
src/fancy_downsa.o: src/fancy_downsa.c src/fancy_image.h
src/fancy_image.o: src/fancy_image.c src/fancy_image.h src/iio.h \
 src/xmalloc.c src/fail.c src/tiff_octaves_rw.c
src/fft.o: src/fft.c src/iio.h src/fail.c src/xmalloc.c src/pickopt.c
src/fftshift.o: src/fftshift.c src/iio.h
src/flambda.o: src/flambda.c src/smapa.h src/fail.c src/xmalloc.c \
 src/random.c src/parsenumbers.c src/colorcoordsf.c src/fancy_image.h \
 src/getpixel.c
src/flowarrows.o: src/flowarrows.c src/iio.h src/fail.c src/xmalloc.c \
 src/drawsegment.c src/getpixel.c src/smapa.h
src/flowinv.o: src/flowinv.c src/iio.h src/fail.c src/xmalloc.c src/bicubic.c \
 src/getpixel.c
src/fontu.o: src/fontu.c src/xmalloc.c src/fail.c src/xfopen.c src/dataconv.c \
 src/iio.h src/pickopt.c
src/getpixel.o: src/getpixel.c
src/ghisto.o: src/ghisto.c src/iio.h src/xmalloc.c src/fail.c src/smapa.h \
 src/pickopt.c
src/gram_schmidt.o: src/gram_schmidt.c
src/grid.o: src/grid.c src/fail.c
src/help_stuff.o: src/help_stuff.c
src/homographies.o: src/homographies.c
src/homwarp.o: src/homwarp.c src/extrapolators.c src/bilinear_interpolation.c \
 src/marching_interpolation.c src/bicubic_gray.c src/spline.c src/iio.h \
 src/xmalloc.c src/fail.c src/parsenumbers.c src/help_stuff.c \
 src/pickopt.c
src/iio.o: src/iio.c
src/iion.o: src/iion.c src/iio.h
src/iion_u16.o: src/iion_u16.c src/iio.h
src/im.o: src/im.c src/all_mains.inc src/ftr/all_mains.inc
src/imflip.o: src/imflip.c src/help_stuff.c src/iio.h
src/imprintf.o: src/imprintf.c src/iio.h src/help_stuff.c
src/lrcat.o: src/lrcat.c src/iio.h src/xmalloc.c src/fail.c src/getpixel.c \
 src/pickopt.c
src/marching_interpolation.o: src/marching_interpolation.c
src/marching_squares.o: src/marching_squares.c
src/mediator.o: src/mediator.c src/iio.h
src/minicg.o: src/minicg.c
src/modes_detector.o: src/modes_detector.c src/smapa.h
src/moistiv_epipolar.o: src/moistiv_epipolar.c src/fail.c
src/morsi.o: src/morsi.c src/iio.h src/help_stuff.c
src/nnint.o: src/nnint.c src/abstract_heap.h src/xmalloc.c src/fail.c \
 src/iio.h src/pickopt.c
src/ntiply.o: src/ntiply.c src/iio.h
src/ok_list.o: src/ok_list.c src/fail.c src/xmalloc.c
src/palette.o: src/palette.c src/fail.c src/xmalloc.c src/xfopen.c \
 src/smapa.h src/iio.h src/pickopt.c src/fontu.c src/dataconv.c \
 src/fonts/xfonts_all.c src/fonts/xfont_10x20.c src/fonts/xfont_12x13ja.c \
 src/fonts/xfont_18x18ja.c src/fonts/xfont_18x18ko.c \
 src/fonts/xfont_4x6.c src/fonts/xfont_5x7.c src/fonts/xfont_5x8.c \
 src/fonts/xfont_6x10.c src/fonts/xfont_6x12.c src/fonts/xfont_6x13.c \
 src/fonts/xfont_6x13B.c src/fonts/xfont_6x13O.c src/fonts/xfont_6x9.c \
 src/fonts/xfont_7x13.c src/fonts/xfont_7x13B.c src/fonts/xfont_7x13O.c \
 src/fonts/xfont_7x14B.c src/fonts/xfont_8x13.c src/fonts/xfont_8x13B.c \
 src/fonts/xfont_8x13O.c src/fonts/xfont_9x15.c src/fonts/xfont_9x15B.c \
 src/fonts/xfont_9x18.c src/fonts/xfont_9x18B.c src/fonts/xfont_canny.c \
 src/fonts/xfont_clR6x12.c src/fonts/xfont_helvR12.c
src/parsenumbers.o: src/parsenumbers.c src/xmalloc.c src/fail.c
src/pickopt.o: src/pickopt.c
src/plambda.o: src/plambda.c src/smapa.h src/fail.c src/xmalloc.c \
 src/random.c src/parsenumbers.c src/colorcoordsf.c src/getpixel.c \
 src/iio.h
src/pview.o: src/pview.c src/iio.h src/fail.c src/xmalloc.c src/xfopen.c \
 src/parsenumbers.c src/drawsegment.c src/pickopt.c src/smapa.h \
 src/random.c
src/qauto.o: src/qauto.c src/iio.h src/help_stuff.c src/pickopt.c
src/qeasy.o: src/qeasy.c src/iio.h
src/random.o: src/random.c
src/ransac.o: src/ransac.c src/fail.c src/xmalloc.c src/xfopen.c src/random.c \
 src/ransac_cases.c src/vvector.h src/homographies.c \
 src/moistiv_epipolar.c src/exterior_algebra.c src/parsenumbers.c
src/ransac_cases.o: src/ransac_cases.c src/vvector.h src/homographies.c \
 src/moistiv_epipolar.c src/fail.c src/exterior_algebra.c
src/redim.o: src/redim.c src/iio.h
src/registration.o: src/registration.c src/iio.h
src/seconds.o: src/seconds.c
src/siftie.o: src/siftie.c src/fail.c src/xmalloc.c src/xfopen.c \
 src/parsenumbers.c src/smapa.h src/ok_list.c src/grid.c src/iio.h
src/siftu.o: src/siftu.c src/siftie.c src/fail.c src/xmalloc.c src/xfopen.c \
 src/parsenumbers.c src/smapa.h src/ok_list.c src/grid.c src/iio.h
src/simpois.o: src/simpois.c src/cleant_cgpois.c src/minicg.c src/smapa.h \
 src/iio.h src/pickopt.c
src/spline.o: src/spline.c
src/srmatch.o: src/srmatch.c src/fail.c src/xmalloc.c src/xfopen.c \
 src/siftie.c src/parsenumbers.c src/smapa.h src/ok_list.c src/grid.c \
 src/iio.h src/ransac.c src/random.c src/ransac_cases.c src/vvector.h \
 src/homographies.c src/moistiv_epipolar.c src/exterior_algebra.c
src/strt.o: src/strt.c src/xmalloc.c src/fail.c src/iio.h src/pickopt.c
src/synflow.o: src/synflow.c src/iio.h src/xmalloc.c src/fail.c \
 src/synflow_core.c src/getpixel.c src/marching_interpolation.c \
 src/vvector.h src/homographies.c src/smapa.h
src/synflow_core.o: src/synflow_core.c src/fail.c src/getpixel.c \
 src/marching_interpolation.c src/vvector.h src/homographies.c \
 src/smapa.h
src/tbcat.o: src/tbcat.c src/iio.h src/smapa.h src/xmalloc.c src/fail.c \
 src/getpixel.c src/pickopt.c
src/tiff_octaves_rw.o: src/tiff_octaves_rw.c
src/tiffu.o: src/tiffu.c
src/upsa.o: src/upsa.c src/iio.h src/fail.c src/marching_squares.c \
 src/marching_interpolation.c src/bicubic.c src/getpixel.c src/pickopt.c
src/veco.o: src/veco.c src/iio.h src/fail.c src/xmalloc.c src/random.c \
 src/help_stuff.c src/pickopt.c
src/vecoh.o: src/vecoh.c src/iio.h src/fail.c src/xmalloc.c src/random.c \
 src/modes_detector.c src/smapa.h src/help_stuff.c src/pickopt.c
src/vecov.o: src/vecov.c src/iio.h src/fail.c src/xmalloc.c src/random.c \
 src/smapa.h src/help_stuff.c src/pickopt.c
src/viewflow.o: src/viewflow.c src/iio.h src/fail.c src/drawsegment.c \
 src/colorcoordsf.c src/marching_squares.c
src/xfopen.o: src/xfopen.c src/fail.c
src/xmalloc.o: src/xmalloc.c src/fail.c
src/ftr/blur.o: src/ftr/blur.c src/ftr/fail.c src/ftr/xmalloc.c src/ftr/smapa.h \
 src/ftr/help_stuff.c src/ftr/parsenumbers.c src/ftr/pickopt.c \
 src/ftr/iio.h
src/ftr/dataconv.o: src/ftr/dataconv.c src/ftr/xmalloc.c src/ftr/fail.c
src/ftr/dosdo.o: src/ftr/dosdo.c src/ftr/iio.h src/ftr/ftr.h src/ftr/fontu.c \
 src/ftr/xmalloc.c src/ftr/fail.c src/ftr/xfopen.c src/ftr/dataconv.c \
 src/ftr/fonts/xfont_9x15.c
src/ftr/egm96.o: src/ftr/egm96.c src/ftr/iio.h
src/ftr/epiview.o: src/ftr/epiview.c src/ftr/iio.h src/ftr/ftr.h src/ftr/fontu.c \
 src/ftr/xmalloc.c src/ftr/fail.c src/ftr/xfopen.c src/ftr/dataconv.c \
 src/ftr/fonts/xfont_9x15.c src/ftr/parsenumbers.c
src/ftr/fail.o: src/ftr/fail.c
src/ftr/fancy_image.o: src/ftr/fancy_image.c src/ftr/fancy_image.h src/ftr/iio.h \
 src/ftr/xmalloc.c src/ftr/fail.c src/ftr/tiff_octaves_rw.c
src/ftr/fontu.o: src/ftr/fontu.c src/ftr/xmalloc.c src/ftr/fail.c \
 src/ftr/xfopen.c src/ftr/dataconv.c src/ftr/iio.h src/ftr/pickopt.c
src/ftr/fpan.o: src/ftr/fpan.c src/ftr/iio.h src/ftr/ftr.h src/ftr/xmalloc.c \
 src/ftr/fail.c src/ftr/smapa.h
src/ftr/fpanflip.o: src/ftr/fpanflip.c src/ftr/iio.h src/ftr/ftr.h \
 src/ftr/xmalloc.c src/ftr/fail.c
src/ftr/fpantiff.o: src/ftr/fpantiff.c src/ftr/tiff_octaves_rw.c src/ftr/ftr.h \
 src/ftr/iio.h src/ftr/pickopt.c
src/ftr/ftr.o: src/ftr/ftr.c src/ftr/ftr_x11.c src/ftr/ftr.h \
 src/ftr/ftr_common_inc.c
src/ftr/ftr_common_inc.o: src/ftr/ftr_common_inc.c
src/ftr/ftr_freeglut.o: src/ftr/ftr_freeglut.c src/ftr/ftr.h \
 src/ftr/ftr_common_inc.c
src/ftr/ftr_mains.o: src/ftr/ftr_mains.c src/ftr/seconds.c src/ftr/ftr.h \
 src/ftr/iio.h
src/ftr/ftr_x11.o: src/ftr/ftr_x11.c src/ftr/ftr.h src/ftr/ftr_common_inc.c
src/ftr/help_stuff.o: src/ftr/help_stuff.c
src/ftr/icrop.o: src/ftr/icrop.c src/ftr/iio.h src/ftr/ftr.h
src/ftr/iio.o: src/ftr/iio.c
src/ftr/marching_interpolation.o: src/ftr/marching_interpolation.c
src/ftr/minisimpois.o: src/ftr/minisimpois.c
src/ftr/parsenumbers.o: src/ftr/parsenumbers.c src/ftr/xmalloc.c src/ftr/fail.c
src/ftr/pickopt.o: src/ftr/pickopt.c
src/ftr/pleview.o: src/ftr/pleview.c src/ftr/tiffu.c src/ftr/ftr.h src/ftr/iio.h
src/ftr/powerkill.o: src/ftr/powerkill.c src/ftr/iio.h src/ftr/ftr.h \
 src/ftr/xmalloc.c src/ftr/fail.c src/ftr/pickopt.c
src/ftr/ppsmooth.o: src/ftr/ppsmooth.c src/ftr/iio.h src/ftr/pickopt.c
src/ftr/random.o: src/ftr/random.c
src/ftr/rpc2.o: src/ftr/rpc2.c src/ftr/xfopen.c src/ftr/fail.c
src/ftr/rpcflip.o: src/ftr/rpcflip.c src/ftr/tiff_octaves_rw.c src/ftr/srtm4o.c \
 src/ftr/rpc2.c src/ftr/xfopen.c src/ftr/fail.c src/ftr/ftr.h \
 src/ftr/iio.h src/ftr/xmalloc.c src/ftr/pickopt.c
src/ftr/seconds.o: src/ftr/seconds.c
src/ftr/srt.o: src/ftr/srt.c src/ftr/ftr.h src/ftr/fontu.c src/ftr/xmalloc.c \
 src/ftr/fail.c src/ftr/xfopen.c src/ftr/dataconv.c \
 src/ftr/fonts/xfont_10x20.c
src/ftr/srtm4o.o: src/ftr/srtm4o.c src/ftr/tiff_octaves_rw.c
src/ftr/strt.o: src/ftr/strt.c src/ftr/xmalloc.c src/ftr/fail.c src/ftr/iio.h \
 src/ftr/pickopt.c
src/ftr/tdip.o: src/ftr/tdip.c src/ftr/iio.h src/ftr/strt.c src/ftr/xmalloc.c \
 src/ftr/fail.c src/ftr/smapa.h src/ftr/random.c src/ftr/pickopt.c
src/ftr/tiff_octaves_rw.o: src/ftr/tiff_octaves_rw.c
src/ftr/tiffu.o: src/ftr/tiffu.c
src/ftr/viho.o: src/ftr/viho.c src/ftr/ftr.h src/ftr/marching_interpolation.c \
 src/ftr/iio.h
src/ftr/vnav.o: src/ftr/vnav.c src/ftr/iio.h src/ftr/tiffu.c src/ftr/ftr.h \
 src/ftr/minisimpois.c src/ftr/blur.c src/ftr/fail.c src/ftr/xmalloc.c \
 src/ftr/smapa.h src/ftr/tdip.c src/ftr/strt.c src/ftr/random.c \
 src/ftr/fontu.c src/ftr/xfopen.c src/ftr/dataconv.c \
 src/ftr/fonts/xfont_9x15.c src/ftr/seconds.c
src/ftr/wifpan.o: src/ftr/wifpan.c src/ftr/iio.h src/ftr/ftr.h src/ftr/xmalloc.c \
 src/ftr/fail.c src/ftr/blur.c src/ftr/smapa.h src/ftr/ppsmooth.c
src/ftr/xfopen.o: src/ftr/xfopen.c src/ftr/fail.c
src/ftr/xmalloc.o: src/ftr/xmalloc.c src/ftr/fail.c
src/misc/abstract_dsf.o: src/misc/abstract_dsf.c
src/misc/aff3d.o: src/misc/aff3d.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/parsenumbers.c src/misc/xmalloc.c
src/misc/alphadots.o: src/misc/alphadots.c src/misc/iio.h src/misc/pickopt.c
src/misc/amle_rec.o: src/misc/amle_rec.c src/misc/smapa.h src/misc/iio.h
src/misc/amle_recsep.o: src/misc/amle_recsep.c src/misc/smapa.h src/misc/iio.h
src/misc/angleplot.o: src/misc/angleplot.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c
src/misc/aronsson11.o: src/misc/aronsson11.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/aronsson43.o: src/misc/aronsson43.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/autholes.o: src/misc/autholes.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/ccproc.c src/misc/abstract_dsf.c
src/misc/awgn.o: src/misc/awgn.c src/misc/smapa.h src/misc/random.c src/misc/iio.h
src/misc/bandhisto.o: src/misc/bandhisto.c src/misc/iio.h
src/misc/bayerdots.o: src/misc/bayerdots.c src/misc/iio.h
src/misc/bayerparts.o: src/misc/bayerparts.c src/misc/iio.h
src/misc/bfollow.o: src/misc/bfollow.c
src/misc/bicubic.o: src/misc/bicubic.c src/misc/getpixel.c
src/misc/bilinear_interpolation.o: src/misc/bilinear_interpolation.c
src/misc/blur.o: src/misc/blur.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/smapa.h src/misc/help_stuff.c src/misc/parsenumbers.c \
 src/misc/pickopt.c src/misc/iio.h
src/misc/bmfm.o: src/misc/bmfm.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/getpixel.c src/misc/iio.h src/misc/parsenumbers.c
src/misc/bmfm_fancier.o: src/misc/bmfm_fancier.c src/misc/xmalloc.c \
 src/misc/fail.c src/misc/getpixel.c src/misc/iio.h \
 src/misc/parsenumbers.c
src/misc/bmfm_fancierw.o: src/misc/bmfm_fancierw.c src/misc/xmalloc.c \
 src/misc/fail.c src/misc/getpixel.c src/misc/iio.h \
 src/misc/parsenumbers.c
src/misc/ccmath_svd.o: src/misc/ccmath_svd.c
src/misc/ccproc.o: src/misc/ccproc.c src/misc/abstract_dsf.c src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/cdeint.o: src/misc/cdeint.c src/misc/iio.h
src/misc/cgi.o: src/misc/cgi.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/bilinear_interpolation.c src/misc/random.c src/misc/smapa.h \
 src/misc/iio.h
src/misc/cglap.o: src/misc/cglap.c src/misc/conjugate_gradient.c \
 src/misc/xmalloc.c src/misc/fail.c src/misc/smapa.h src/misc/iio.h
src/misc/cgpois.o: src/misc/cgpois.c src/misc/conjugate_gradient.c \
 src/misc/xmalloc.c src/misc/fail.c src/misc/smapa.h src/misc/iio.h
src/misc/cgpois_rec.o: src/misc/cgpois_rec.c src/misc/conjugate_gradient.c \
 src/misc/xmalloc.c src/misc/fail.c src/misc/smapa.h src/misc/iio.h
src/misc/chisto.o: src/misc/chisto.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/smapa.h
src/misc/cldmask.o: src/misc/cldmask.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/xfopen.c src/misc/parsenumbers.c src/misc/drawsegment.c \
 src/misc/iio.h src/misc/pickopt.c
src/misc/cline.o: src/misc/cline.c src/misc/fail.c src/misc/bicubic.c \
 src/misc/getpixel.c src/misc/bilinear_interpolation.c src/misc/iio.h \
 src/misc/smapa.h
src/misc/closeup.o: src/misc/closeup.c src/misc/iio.h src/misc/marching_squares.c \
 src/misc/marching_interpolation.c src/misc/bicubic.c src/misc/getpixel.c
src/misc/colorcoords.o: src/misc/colorcoords.c
src/misc/colorflow.o: src/misc/colorflow.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/colorcoords.c src/misc/drawsegment.c \
 src/misc/marching_squares.c src/misc/iio.h
src/misc/colormesh.o: src/misc/colormesh.c src/misc/iio.h src/misc/fail.c \
 src/misc/rpc.c src/misc/xfopen.c src/misc/smapa.h
src/misc/colormeshh.o: src/misc/colormeshh.c src/misc/iio.h src/misc/fail.c \
 src/misc/rpc.c src/misc/xfopen.c src/misc/smapa.h
src/misc/component.o: src/misc/component.c src/misc/iio.h
src/misc/conjugate_gradient.o: src/misc/conjugate_gradient.c src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/conjugate_gradient_linverif.o: src/misc/conjugate_gradient_linverif.c \
 src/misc/random.c
src/misc/contours.o: src/misc/contours.c src/misc/marching_squares.c \
 src/misc/marching_interpolation.c src/misc/iio.h
src/misc/convolution.o: src/misc/convolution.c
src/misc/crop_rect.o: src/misc/crop_rect.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/iio.h
src/misc/crop_size.o: src/misc/crop_size.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/iio.h
src/misc/croparound.o: src/misc/croparound.c src/misc/getpixel.c src/misc/iio.h \
 src/misc/xmalloc.c src/misc/fail.c
src/misc/crosses.o: src/misc/crosses.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/iio.h
src/misc/cutrecombine.o: src/misc/cutrecombine.c src/misc/iio.h src/misc/fail.c
src/misc/dataconv.o: src/misc/dataconv.c src/misc/xmalloc.c src/misc/fail.c
src/misc/deframe.o: src/misc/deframe.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/getpixel.c src/misc/homographies.c src/misc/bicubic.c \
 src/misc/iio.h
src/misc/disp_to_corr.o: src/misc/disp_to_corr.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c
src/misc/distance.o: src/misc/distance.c src/misc/abstract_heap.h \
 src/misc/xmalloc.c src/misc/fail.c src/misc/iio.h src/misc/xfopen.c \
 src/misc/parsenumbers.c src/misc/smapa.h
src/misc/ditheru.o: src/misc/ditheru.c src/misc/iio.h src/misc/pickopt.c
src/misc/drawroads.o: src/misc/drawroads.c src/misc/drawsegment.c src/misc/iio.h \
 src/misc/fail.c src/misc/parsenumbers.c src/misc/xmalloc.c \
 src/misc/pickopt.c
src/misc/drawsegment.o: src/misc/drawsegment.c
src/misc/dveco.o: src/misc/dveco.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/random.c src/misc/pickopt.c
src/misc/dvecov_lm.o: src/misc/dvecov_lm.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/pickopt.c
src/misc/elap.o: src/misc/elap.c src/misc/iio.h
src/misc/elap2.o: src/misc/elap2.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/distance.c src/misc/abstract_heap.h src/misc/iio.h
src/misc/elap3.o: src/misc/elap3.c src/misc/iio.h
src/misc/elap_rec.o: src/misc/elap_rec.c src/misc/iio.h
src/misc/elap_recsep.o: src/misc/elap_recsep.c src/misc/smapa.h src/misc/iio.h
src/misc/elevate_matches.o: src/misc/elevate_matches.c src/misc/iio.h \
 src/misc/xfopen.c src/misc/fail.c
src/misc/elevate_matcheshh.o: src/misc/elevate_matcheshh.c src/misc/iio.h \
 src/misc/xfopen.c src/misc/fail.c src/misc/parsenumbers.c \
 src/misc/xmalloc.c
src/misc/epilist.o: src/misc/epilist.c
src/misc/error.o: src/misc/error.c
src/misc/fabius.o: src/misc/fabius.c
src/misc/fail.o: src/misc/fail.c
src/misc/fancy_evals.o: src/misc/fancy_evals.c src/misc/fancy_image.h
src/misc/fancy_image.o: src/misc/fancy_image.c src/misc/fancy_image.h \
 src/misc/iio.h src/misc/xmalloc.c src/misc/fail.c \
 src/misc/tiff_octaves_rw.c
src/misc/fancy_zoomout.o: src/misc/fancy_zoomout.c src/misc/fancy_image.h
src/misc/faxpb.o: src/misc/faxpb.c src/misc/iio.h
src/misc/faxpby.o: src/misc/faxpby.c src/misc/iio.h
src/misc/fft.o: src/misc/fft.c src/misc/iio.h src/misc/fail.c src/misc/xmalloc.c \
 src/misc/pickopt.c
src/misc/fftper.o: src/misc/fftper.c src/misc/iio.h
src/misc/fftwisdom.o: src/misc/fftwisdom.c
src/misc/fill_bill.o: src/misc/fill_bill.c src/misc/iio.h
src/misc/fill_rect.o: src/misc/fill_rect.c src/misc/iio.h
src/misc/fillcorners.o: src/misc/fillcorners.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/flow_ms.o: src/misc/flow_ms.c src/misc/iio.h src/misc/smapa.h \
 src/misc/fail.c src/misc/xmalloc.c
src/misc/flowback.o: src/misc/flowback.c src/misc/iio.h
src/misc/flowdiv.o: src/misc/flowdiv.c src/misc/iio.h src/misc/fragments.c \
 src/misc/getpixel.c
src/misc/flowgrad.o: src/misc/flowgrad.c src/misc/iio.h src/misc/fragments.c \
 src/misc/getpixel.c
src/misc/flowh.o: src/misc/flowh.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/pickopt.c
src/misc/flowjac.o: src/misc/flowjac.c src/misc/iio.h src/misc/fragments.c \
 src/misc/getpixel.c
src/misc/flownop.o: src/misc/flownop.c src/misc/iio.h src/misc/fragments.c \
 src/misc/getpixel.c
src/misc/fmsrA.o: src/misc/fmsrA.c src/misc/svd.c
src/misc/fnorm.o: src/misc/fnorm.c src/misc/iio.h
src/misc/font_6x12.o: src/misc/font_6x12.c
src/misc/fpgraph.o: src/misc/fpgraph.c src/misc/iio.h src/misc/pickopt.c
src/misc/fragments.o: src/misc/fragments.c
src/misc/frakes_monaco_smith.o: src/misc/frakes_monaco_smith.c src/misc/fail.c \
 src/misc/xmalloc.c src/misc/iio.h
src/misc/frustumize.o: src/misc/frustumize.c src/misc/iio.h \
 src/misc/parsenumbers.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/pickopt.c src/misc/smapa.h
src/misc/gblur.o: src/misc/gblur.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/vvector.h
src/misc/gblur_core.o: src/misc/gblur_core.c src/misc/iio.h
src/misc/genk.o: src/misc/genk.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/random.c src/misc/smapa.h src/misc/iio.h
src/misc/getopt_with_pickopt.o: src/misc/getopt_with_pickopt.c
src/misc/getpixel.o: src/misc/getpixel.c
src/misc/gharrows.o: src/misc/gharrows.c src/misc/iio.h src/misc/smapa.h
src/misc/ghough.o: src/misc/ghough.c src/misc/iio.h
src/misc/ghough2.o: src/misc/ghough2.c src/misc/pickopt.c src/misc/iio.h
src/misc/graysing.o: src/misc/graysing.c src/misc/iio.h
src/misc/grid.o: src/misc/grid.c src/misc/fail.c
src/misc/harris.o: src/misc/harris.c src/misc/iio.h src/misc/fragments.c \
 src/misc/getpixel.c
src/misc/help_stuff.o: src/misc/help_stuff.c
src/misc/histeq8.o: src/misc/histeq8.c src/misc/iio.h
src/misc/histomodev.o: src/misc/histomodev.c src/misc/iio.h src/misc/pickopt.c
src/misc/homaff.o: src/misc/homaff.c
src/misc/homdots.o: src/misc/homdots.c src/misc/iio.h
src/misc/homfilt.o: src/misc/homfilt.c src/misc/parsenumbers.c src/misc/xmalloc.c \
 src/misc/fail.c src/misc/vvector.h src/misc/smapa.h
src/misc/homographies.o: src/misc/homographies.c
src/misc/houghs.o: src/misc/houghs.c src/misc/iio.h
src/misc/hrezoom.o: src/misc/hrezoom.c src/misc/vvector.h
src/misc/hs.o: src/misc/hs.c src/misc/iio.h
src/misc/hs_core.o: src/misc/hs_core.c
src/misc/hs_new.o: src/misc/hs_new.c src/misc/iio.h
src/misc/huffman.o: src/misc/huffman.c src/misc/abstract_heap.h
src/misc/hview.o: src/misc/hview.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/iio.h
src/misc/ihough2.o: src/misc/ihough2.c src/misc/iio.h src/misc/pickopt.c
src/misc/iio.o: src/misc/iio.c
src/misc/ijmesh.o: src/misc/ijmesh.c src/misc/iio.h
src/misc/imdim.o: src/misc/imdim.c src/misc/iio.h
src/misc/imgerr.o: src/misc/imgerr.c src/misc/fail.c src/misc/iio.h \
 src/misc/smapa.h src/misc/xmalloc.c
src/misc/imgstats.o: src/misc/imgstats.c src/misc/iio.h
src/misc/iminfo.o: src/misc/iminfo.c src/misc/iio.h src/misc/statistics.c \
 src/misc/fail.c src/misc/xmalloc.c
src/misc/imspread.o: src/misc/imspread.c src/misc/iio.h
src/misc/inppairs.o: src/misc/inppairs.c src/misc/iio.h
src/misc/intimg.o: src/misc/intimg.c src/misc/iio.h
src/misc/ipol_datum.o: src/misc/ipol_datum.c
src/misc/ipol_watermark.o: src/misc/ipol_watermark.c src/misc/iio.h
src/misc/ising.o: src/misc/ising.c src/misc/random.c
src/misc/isingroot.o: src/misc/isingroot.c src/misc/random.c src/misc/iio.h \
 src/misc/pickopt.c
src/misc/isoricci.o: src/misc/isoricci.c src/misc/iio.h
src/misc/lapbediag.o: src/misc/lapbediag.c src/misc/smapa.h src/misc/iio.h
src/misc/lapbediag_sep.o: src/misc/lapbediag_sep.c src/misc/smapa.h src/misc/iio.h
src/misc/lapcolo.o: src/misc/lapcolo.c src/misc/lapbediag.c src/misc/smapa.h \
 src/misc/iio.h
src/misc/lgblur.o: src/misc/lgblur.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c
src/misc/lgblur2.o: src/misc/lgblur2.c src/misc/gblur.c src/misc/iio.h \
 src/misc/fail.c src/misc/xmalloc.c src/misc/vvector.h
src/misc/lgblur3.o: src/misc/lgblur3.c src/misc/gblur.c src/misc/iio.h \
 src/misc/fail.c src/misc/xmalloc.c src/misc/vvector.h src/misc/bicubic.c \
 src/misc/getpixel.c src/misc/smapa.h
src/misc/lic.o: src/misc/lic.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/random.c src/misc/bilinear_interpolation.c src/misc/smapa.h \
 src/misc/iio.h
src/misc/lk.o: src/misc/lk.c src/misc/iio.h src/misc/svd.c src/misc/vvector.h \
 src/misc/smapa.h
src/misc/lk_omp.o: src/misc/lk_omp.c src/misc/iio.h src/misc/svd.c \
 src/misc/vvector.h
src/misc/lkgen.o: src/misc/lkgen.c
src/misc/lure.o: src/misc/lure.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/getpixel.c src/misc/blur.c src/misc/smapa.h src/misc/iio.h
src/misc/lures.o: src/misc/lures.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/getpixel.c src/misc/blur.c src/misc/smapa.h src/misc/iio.h
src/misc/maptp.o: src/misc/maptp.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/parsenumbers.c
src/misc/marchi_clean.o: src/misc/marchi_clean.c
src/misc/marching_follow.o: src/misc/marching_follow.c
src/misc/marching_interpolation.o: src/misc/marching_interpolation.c
src/misc/marching_squares.o: src/misc/marching_squares.c
src/misc/mediatorw.o: src/misc/mediatorw.c src/misc/iio.h
src/misc/metatiler.o: src/misc/metatiler.c src/misc/iio.h
src/misc/mima.o: src/misc/mima.c src/misc/iio.h
src/misc/minimize.o: src/misc/minimize.c src/misc/fail.c src/misc/minimize_gsl.c \
 src/misc/smapa.h
src/misc/minimize_ccmath.o: src/misc/minimize_ccmath.c src/misc/optmiz.c \
 src/misc/smapa.h
src/misc/minimize_gsl.o: src/misc/minimize_gsl.c src/misc/smapa.h
src/misc/mnehs.o: src/misc/mnehs.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/getpixel.c src/misc/bicubic.c src/misc/rpc.c \
 src/misc/xfopen.c src/misc/smapa.h src/misc/parsenumbers.c \
 src/misc/pickopt.c
src/misc/mnehs_ms.o: src/misc/mnehs_ms.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/getpixel.c src/misc/bicubic.c src/misc/smapa.h \
 src/misc/parsenumbers.c src/misc/pickopt.c
src/misc/monodisp.o: src/misc/monodisp.c src/misc/iio.h
src/misc/morsi_demo.o: src/misc/morsi_demo.c src/misc/iio.h
src/misc/morsi_demo2.o: src/misc/morsi_demo2.c src/misc/iio.h
src/misc/ofc.o: src/misc/ofc.c src/misc/iio.h src/misc/fragments.c \
 src/misc/vvector.h
src/misc/ok_list.o: src/misc/ok_list.c src/misc/fail.c src/misc/xmalloc.c
src/misc/optmiz.o: src/misc/optmiz.c
src/misc/overflow.o: src/misc/overflow.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/xfopen.c
src/misc/overpoints.o: src/misc/overpoints.c src/misc/iio.h src/misc/fail.c \
 src/misc/xfopen.c src/misc/xmalloc.c src/misc/smapa.h
src/misc/p.o: src/misc/p.c
src/misc/pairhom.o: src/misc/pairhom.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/parsenumbers.c src/misc/pickopt.c
src/misc/pairsinp.o: src/misc/pairsinp.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/parsenumbers.c
src/misc/pamle_rec.o: src/misc/pamle_rec.c src/misc/smapa.h src/misc/iio.h
src/misc/paraflow.o: src/misc/paraflow.c src/misc/iio.h src/misc/fail.c \
 src/misc/statistics.c src/misc/xmalloc.c src/misc/synflow_core.c \
 src/misc/getpixel.c src/misc/marching_interpolation.c src/misc/vvector.h \
 src/misc/homographies.c src/misc/smapa.h
src/misc/parsenumbers.o: src/misc/parsenumbers.c src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/periodize.o: src/misc/periodize.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/perms.o: src/misc/perms.c
src/misc/pickopt.o: src/misc/pickopt.c
src/misc/plyflatten.o: src/misc/plyflatten.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/smapa.h src/misc/iio.h
src/misc/plyroads.o: src/misc/plyroads.c src/misc/fail.c src/misc/xfopen.c \
 src/misc/parsenumbers.c src/misc/xmalloc.c
src/misc/plyroads_mini.o: src/misc/plyroads_mini.c src/misc/parsenumbers.c \
 src/misc/xmalloc.c src/misc/fail.c
src/misc/pmba.o: src/misc/pmba.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/parsenumbers.c src/misc/xmalloc.c src/misc/pickopt.c
src/misc/pmba2.o: src/misc/pmba2.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/parsenumbers.c src/misc/xmalloc.c src/misc/pickopt.c
src/misc/poisson_rec.o: src/misc/poisson_rec.c src/misc/iio.h
src/misc/polygonify.o: src/misc/polygonify.c src/misc/ccproc.c \
 src/misc/abstract_dsf.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/smapa.h src/misc/iio.h
src/misc/pomini.o: src/misc/pomini.c
src/misc/posmax.o: src/misc/posmax.c src/misc/iio.h
src/misc/ppsboundary.o: src/misc/ppsboundary.c src/misc/iio.h
src/misc/ppsmooth.o: src/misc/ppsmooth.c src/misc/iio.h src/misc/pickopt.c
src/misc/printmask.o: src/misc/printmask.c src/misc/iio.h
src/misc/q.o: src/misc/q.c
src/misc/raddots.o: src/misc/raddots.c src/misc/iio.h
src/misc/radphar.o: src/misc/radphar.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/random.o: src/misc/random.c
src/misc/ranrecs.o: src/misc/ranrecs.c src/misc/iio.h src/misc/fragments.c
src/misc/really_simplest_inpainting.o: src/misc/really_simplest_inpainting.c \
 src/misc/iio.h
src/misc/reallysmall_svd.o: src/misc/reallysmall_svd.c
src/misc/reallysmall_svd2.o: src/misc/reallysmall_svd2.c
src/misc/reallysmall_svd3.o: src/misc/reallysmall_svd3.c
src/misc/rednoise.o: src/misc/rednoise.c src/misc/iio.h
src/misc/remove_small_cc.o: src/misc/remove_small_cc.c src/misc/abstract_dsf.c \
 src/misc/iio.h src/misc/pickopt.c
src/misc/replicate.o: src/misc/replicate.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/rgbcube.o: src/misc/rgbcube.c src/misc/iio.h
src/misc/rgfield.o: src/misc/rgfield.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/smapa.h src/misc/gblur.c src/misc/vvector.h
src/misc/rgfields.o: src/misc/rgfields.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/smapa.h src/misc/gblur.c src/misc/vvector.h
src/misc/rgfieldst.o: src/misc/rgfieldst.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/smapa.h src/misc/gblur.c src/misc/vvector.h
src/misc/rip.o: src/misc/rip.c
src/misc/rpc.o: src/misc/rpc.c src/misc/xfopen.c src/misc/fail.c src/misc/smapa.h
src/misc/rpc_angpair.o: src/misc/rpc_angpair.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/fail.c src/misc/smapa.h src/misc/xmalloc.c src/misc/iio.h
src/misc/rpc_curvpair.o: src/misc/rpc_curvpair.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/fail.c src/misc/smapa.h src/misc/xmalloc.c src/misc/iio.h
src/misc/rpc_epicyl.o: src/misc/rpc_epicyl.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/fail.c src/misc/smapa.h
src/misc/rpc_errfilt.o: src/misc/rpc_errfilt.c src/misc/parsenumbers.c \
 src/misc/xmalloc.c src/misc/fail.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/smapa.h src/misc/pickopt.c
src/misc/rpc_errpair.o: src/misc/rpc_errpair.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/fail.c src/misc/smapa.h src/misc/xmalloc.c src/misc/iio.h
src/misc/rpc_errsingle.o: src/misc/rpc_errsingle.c src/misc/rpc.c \
 src/misc/xfopen.c src/misc/fail.c src/misc/smapa.h src/misc/xmalloc.c \
 src/misc/iio.h
src/misc/rpc_eval.o: src/misc/rpc_eval.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/fail.c src/misc/smapa.h
src/misc/rpc_mnehs.o: src/misc/rpc_mnehs.c src/misc/iio.h src/misc/getpixel.c \
 src/misc/bicubic.c src/misc/rpc.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/smapa.h src/misc/tiffu.c src/misc/xmalloc.c
src/misc/rpc_pm.o: src/misc/rpc_pm.c src/misc/iio.h src/misc/getpixel.c \
 src/misc/bicubic.c src/misc/rpc.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/smapa.h src/misc/tiffu.c src/misc/xmalloc.c
src/misc/rpc_pmn.o: src/misc/rpc_pmn.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/fail.c src/misc/smapa.h src/misc/tiffu.c src/misc/xmalloc.c \
 src/misc/iio.h
src/misc/rpc_warpab.o: src/misc/rpc_warpab.c src/misc/getpixel.c \
 src/misc/bicubic.c src/misc/rpc.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/smapa.h src/misc/iio.h src/misc/xmalloc.c
src/misc/rpc_warpabt.o: src/misc/rpc_warpabt.c src/misc/getpixel.c \
 src/misc/bicubic.c src/misc/rpc.c src/misc/xfopen.c src/misc/fail.c \
 src/misc/smapa.h src/misc/tiffu.c src/misc/iio.h src/misc/xmalloc.c
src/misc/rpcflow.o: src/misc/rpcflow.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c src/misc/rpc.c src/misc/xfopen.c src/misc/smapa.h
src/misc/rpchfilt.o: src/misc/rpchfilt.c src/misc/parsenumbers.c \
 src/misc/xmalloc.c src/misc/fail.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/smapa.h
src/misc/rpcparcheck.o: src/misc/rpcparcheck.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/parsenumbers.c src/misc/rpc.c src/misc/xfopen.c \
 src/misc/smapa.h
src/misc/satproj.o: src/misc/satproj.c src/misc/fail.c src/misc/getpixel.c \
 src/misc/iio.h src/misc/xmalloc.c src/misc/parsenumbers.c \
 src/misc/pickopt.c
src/misc/scheme_plap.o: src/misc/scheme_plap.c src/misc/fail.c src/misc/iio.h \
 src/misc/xmalloc.c
src/misc/sdistance.o: src/misc/sdistance.c src/misc/abstract_heap.h \
 src/misc/xmalloc.c src/misc/fail.c src/misc/iio.h src/misc/pickopt.c
src/misc/segfilter.o: src/misc/segfilter.c
src/misc/setdim.o: src/misc/setdim.c src/misc/iio.h
src/misc/sfblur.o: src/misc/sfblur.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/iio.h
src/misc/shadowcast.o: src/misc/shadowcast.c src/misc/iio.h src/misc/pickopt.c
src/misc/sheartilter.o: src/misc/sheartilter.c
src/misc/shuntingyard.o: src/misc/shuntingyard.c src/misc/fail.c \
 src/misc/xmalloc.c
src/misc/simplest_inpainting.o: src/misc/simplest_inpainting.c src/misc/iio.h
src/misc/sort.o: src/misc/sort.c src/misc/iio.h
src/misc/sphereheights.o: src/misc/sphereheights.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c
src/misc/spline_double.o: src/misc/spline_double.c
src/misc/starfield.o: src/misc/starfield.c src/misc/random.c src/misc/iio.h \
 src/misc/xmalloc.c src/misc/fail.c
src/misc/startracker.o: src/misc/startracker.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/getpixel.c src/misc/ok_list.c src/misc/grid.c src/misc/iio.h
src/misc/statistics.o: src/misc/statistics.c src/misc/fail.c src/misc/xmalloc.c
src/misc/sterint.o: src/misc/sterint.c src/misc/iio.h src/misc/xmalloc.c \
 src/misc/fail.c
src/misc/strt.o: src/misc/strt.c src/misc/xmalloc.c src/misc/fail.c src/misc/iio.h \
 src/misc/pickopt.c
src/misc/sum.o: src/misc/sum.c src/misc/iio.h
src/misc/svd.o: src/misc/svd.c
src/misc/svd_example.o: src/misc/svd_example.c
src/misc/svd_smallexample.o: src/misc/svd_smallexample.c \
 src/misc/reallysmall_svd.c
src/misc/svdmini.o: src/misc/svdmini.c
src/misc/synflow_core.o: src/misc/synflow_core.c src/misc/fail.c \
 src/misc/getpixel.c src/misc/marching_interpolation.c src/misc/vvector.h \
 src/misc/homographies.c src/misc/smapa.h
src/misc/tcregistration.o: src/misc/tcregistration.c src/misc/iio.h
src/misc/tdip.o: src/misc/tdip.c src/misc/iio.h src/misc/strt.c src/misc/xmalloc.c \
 src/misc/fail.c src/misc/smapa.h src/misc/random.c src/misc/pickopt.c
src/misc/testgblur.o: src/misc/testgblur.c src/misc/gblur.c src/misc/iio.h \
 src/misc/fail.c src/misc/xmalloc.c src/misc/vvector.h
src/misc/thinpa_recsep.o: src/misc/thinpa_recsep.c src/misc/smapa.h src/misc/iio.h
src/misc/tiff_octaves.o: src/misc/tiff_octaves.c
src/misc/tiff_octaves_notest.o: src/misc/tiff_octaves_notest.c
src/misc/tiff_octaves_old.o: src/misc/tiff_octaves_old.c
src/misc/tiff_octaves_rw.o: src/misc/tiff_octaves_rw.c
src/misc/tiffu.o: src/misc/tiffu.c
src/misc/tiloct.o: src/misc/tiloct.c
src/misc/tilt_and_shear.o: src/misc/tilt_and_shear.c src/misc/bicubic.c \
 src/misc/getpixel.c src/misc/fail.c src/misc/xmalloc.c src/misc/iio.h
src/misc/tiny_lure.o: src/misc/tiny_lure.c src/misc/getpixel.c src/misc/blur.c \
 src/misc/fail.c src/misc/xmalloc.c src/misc/smapa.h src/misc/iio.h
src/misc/tregistration.o: src/misc/tregistration.c src/misc/iio.h
src/misc/try_fancy.o: src/misc/try_fancy.c src/misc/fancy_image.h
src/misc/tvint.o: src/misc/tvint.c src/misc/smapa.h src/misc/iio.h \
 src/misc/pickopt.c
src/misc/unalpha.o: src/misc/unalpha.c src/misc/iio.h
src/misc/uncrop.o: src/misc/uncrop.c src/misc/iio.h
src/misc/unshadow.o: src/misc/unshadow.c src/misc/iio.h src/misc/pickopt.c
src/misc/utm.o: src/misc/utm.c
src/misc/vecov_lm.o: src/misc/vecov_lm.c src/misc/iio.h src/misc/fail.c \
 src/misc/xmalloc.c src/misc/pickopt.c
src/misc/vecstack.o: src/misc/vecstack.c src/misc/iio.h
src/misc/watermark.o: src/misc/watermark.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/xfopen.c src/misc/font_6x12.c src/misc/iio.h
src/misc/wgs84.o: src/misc/wgs84.c
src/misc/xfopen.o: src/misc/xfopen.c src/misc/fail.c
src/misc/xmalloc.o: src/misc/xmalloc.c src/misc/fail.c
src/misc/xmalloc_stats.o: src/misc/xmalloc_stats.c
src/misc/zeropad.o: src/misc/zeropad.c src/misc/xmalloc.c src/misc/fail.c \
 src/misc/iio.h
src/misc/zoombil.o: src/misc/zoombil.c src/misc/fail.c src/misc/xmalloc.c \
 src/misc/bilinear_interpolation.c src/misc/smapa.h src/misc/iio.h
