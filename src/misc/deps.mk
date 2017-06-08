aff3d.o: aff3d.c xfopen.c fail.c parsenumbers.c xmalloc.c
alphadots.o: alphadots.c iio.h pickopt.c
amle_rec.o: amle_rec.c smapa.h iio.h
amle_recsep.o: amle_recsep.c smapa.h iio.h
angleplot.o: angleplot.c iio.h fail.c xmalloc.c
apm.o: apm.c iio.h smapa.h
aronsson11.o: aronsson11.c iio.h xmalloc.c fail.c
aronsson43.o: aronsson43.c iio.h xmalloc.c fail.c
autholes.o: autholes.c iio.h fail.c xmalloc.c ccproc.c abstract_dsf.c
awgn.o: awgn.c smapa.h random.c iio.h
bandhisto.o: bandhisto.c iio.h
basify.o: basify.c
bayerdots.o: bayerdots.c iio.h
bayerparts.o: bayerparts.c iio.h
blur.o: blur.c fail.c xmalloc.c smapa.h iio.h parsenumbers.c
bmfm.o: bmfm.c xmalloc.c fail.c getpixel.c
bmfm_fancier.o: bmfm_fancier.c xmalloc.c fail.c getpixel.c
bmfm_fancierw.o: bmfm_fancierw.c xmalloc.c fail.c getpixel.c
ccfilt.o: ccfilt.c xmalloc.c fail.c iio.h
ccproc.o: ccproc.c abstract_dsf.c xmalloc.c fail.c
cdeint.o: cdeint.c iio.h
centroid.o: centroid.c bicubic.c getpixel.c iio.h
cgi.o: cgi.c fail.c xmalloc.c bilinear_interpolation.c random.c smapa.h \
 iio.h
cglap.o: cglap.c conjugate_gradient.c xmalloc.c fail.c smapa.h iio.h
cglapsep.o: cglapsep.c conjugate_gradient.c xmalloc.c fail.c smapa.h \
 iio.h
cgpois.o: cgpois.c conjugate_gradient.c xmalloc.c fail.c smapa.h iio.h
cgpois_rec.o: cgpois_rec.c conjugate_gradient.c xmalloc.c fail.c smapa.h \
 iio.h
chan.o: chan.c iio.h
chisto.o: chisto.c iio.h xmalloc.c fail.c smapa.h
cldmask.o: cldmask.c fail.c xmalloc.c xfopen.c parsenumbers.c \
 drawsegment.c iio.h pickopt.c
cline.o: cline.c fail.c bicubic.c getpixel.c bilinear_interpolation.c \
 iio.h smapa.h
closeup.o: closeup.c iio.h marching_squares.c marching_interpolation.c \
 bicubic.c getpixel.c
colorflow.o: colorflow.c fail.c xmalloc.c colorcoords.c drawsegment.c \
 marching_squares.c iio.h
colormesh.o: colormesh.c iio.h fail.c rpc.c xfopen.c smapa.h
colormeshh.o: colormeshh.c iio.h fail.c rpc.c xfopen.c smapa.h
component.o: component.c iio.h
contours.o: contours.c marching_squares.c marching_interpolation.c iio.h
convolution.o: convolution.c
crop_rect.o: crop_rect.c fail.c xmalloc.c iio.h
crop_size.o: crop_size.c fail.c xmalloc.c iio.h
croparound.o: croparound.c getpixel.c iio.h xmalloc.c fail.c
crosses.o: crosses.c fail.c xmalloc.c iio.h
cutrecombine.o: cutrecombine.c iio.h fail.c
darboux.o: darboux.c xmalloc.c fail.c iio.h
deframe.o: deframe.c fail.c xmalloc.c getpixel.c homographies.c bicubic.c \
 iio.h
disp_to_corr.o: disp_to_corr.c iio.h fail.c xmalloc.c
distance.o: distance.c abstract_heap.h xmalloc.c fail.c iio.h xfopen.c \
 parsenumbers.c smapa.h
ditheru.o: ditheru.c iio.h pickopt.c
drawroads.o: drawroads.c drawsegment.c iio.h fail.c parsenumbers.c \
 xmalloc.c pickopt.c
dveco.o: dveco.c iio.h fail.c xmalloc.c random.c pickopt.c
dvecov_lm.o: dvecov_lm.c iio.h fail.c xmalloc.c pickopt.c
elap.o: elap.c iio.h
elap2.o: elap2.c iio.h
elap3.o: elap3.c iio.h
elap_rec.o: elap_rec.c iio.h
elap_recsep.o: elap_recsep.c smapa.h iio.h
elevate_matches.o: elevate_matches.c iio.h xfopen.c fail.c
elevate_matcheshh.o: elevate_matcheshh.c iio.h xfopen.c fail.c \
 parsenumbers.c xmalloc.c
fabius.o: fabius.c
fancy_evals.o: fancy_evals.c fancy_image.c fancy_image.h iio.h xmalloc.c \
 fail.c tiff_octaves_rw.c
fancy_image.o: fancy_image.c fancy_image.h iio.h xmalloc.c fail.c \
 tiff_octaves_rw.c
fancy_zoomout.o: fancy_zoomout.c fancy_image.h
faxpb.o: faxpb.c iio.h
faxpby.o: faxpby.c iio.h
fft.o: fft.c iio.h fail.c xmalloc.c
fftper.o: fftper.c iio.h
fill_bill.o: fill_bill.c iio.h
fill_rect.o: fill_rect.c iio.h
fillcorners.o: fillcorners.c iio.h xmalloc.c fail.c
flow_ms.o: flow_ms.c iio.h smapa.h fail.c xmalloc.c
flowback.o: flowback.c iio.h
flowdiv.o: flowdiv.c iio.h fragments.c getpixel.c
flowgrad.o: flowgrad.c iio.h fragments.c getpixel.c
flowh.o: flowh.c iio.h fail.c xmalloc.c
flowjac.o: flowjac.c iio.h fragments.c getpixel.c
flownop.o: flownop.c iio.h fragments.c getpixel.c
fmsrA.o: fmsrA.c svd.c
fnorm.o: fnorm.c iio.h
fontu2.o: fontu2.c xmalloc.c fail.c xfopen.c dataconv.c iio.h pickopt.c
fpgraph.o: fpgraph.c iio.h pickopt.c
frakes_monaco_smith.o: frakes_monaco_smith.c fail.c xmalloc.c iio.h
frustumize.o: frustumize.c iio.h parsenumbers.c xmalloc.c fail.c \
 pickopt.c smapa.h
gblur.o: gblur.c iio.h fail.c xmalloc.c vvector.h
gblur_core.o: gblur_core.c iio.h
gblurs.o: gblurs.c iio.h
genk.o: genk.c fail.c xmalloc.c random.c smapa.h iio.h
gharrows.o: gharrows.c iio.h smapa.h
ghough.o: ghough.c iio.h
ghough2.o: ghough2.c pickopt.c iio.h
graysing.o: graysing.c iio.h
harris.o: harris.c iio.h fragments.c getpixel.c
histeq8.o: histeq8.c iio.h
histomodev.o: histomodev.c iio.h pickopt.c
homdots.o: homdots.c iio.h
homfilt.o: homfilt.c parsenumbers.c xmalloc.c fail.c vvector.h smapa.h
houghs.o: houghs.c iio.h
hrezoom.o: hrezoom.c vvector.h
hs.o: hs.c iio.h
hs_init.o: hs_init.c iio.h
hs_omp.o: hs_omp.c iio.h
huffman.o: huffman.c abstract_heap.h
hview.o: hview.c xmalloc.c fail.c iio.h
ihough2.o: ihough2.c iio.h pickopt.c
ijmesh.o: ijmesh.c iio.h
imdim.o: imdim.c iio.h
imgdist.o: imgdist.c iio.h
imgerr.o: imgerr.c fail.c iio.h smapa.h xmalloc.c
imgstats.o: imgstats.c iio.h
iminfo.o: iminfo.c iio.h statistics.c fail.c xmalloc.c
imspread.o: imspread.c iio.h
inppairs.o: inppairs.c iio.h
intimg.o: intimg.c iio.h
ipol_watermark.o: ipol_watermark.c iio.h
isingroot.o: isingroot.c random.c iio.h pickopt.c
isoricci.o: isoricci.c iio.h
lapbediag.o: lapbediag.c smapa.h iio.h
lapbediag_sep.o: lapbediag_sep.c smapa.h iio.h
lapcolo.o: lapcolo.c lapbediag.c smapa.h iio.h
lgblur.o: lgblur.c iio.h fail.c xmalloc.c
lgblur2.o: lgblur2.c gblur.c iio.h fail.c xmalloc.c vvector.h
lgblur3.o: lgblur3.c gblur.c iio.h fail.c xmalloc.c vvector.h bicubic.c \
 getpixel.c smapa.h
lic.o: lic.c fail.c xmalloc.c random.c bilinear_interpolation.c smapa.h \
 iio.h
lk.o: lk.c iio.h svd.c vvector.h smapa.h
lk_omp.o: lk_omp.c iio.h svd.c vvector.h
lka.o: lka.c iio.h svd.c vvector.h
lure.o: lure.c xmalloc.c fail.c getpixel.c blur.c smapa.h iio.h
lures.o: lures.c xmalloc.c fail.c getpixel.c blur.c smapa.h iio.h
maptp.o: maptp.c fail.c xmalloc.c parsenumbers.c
marching_squares.o: marching_squares.c
metatiler.o: metatiler.c iio.h
mima.o: mima.c iio.h
mini_fhtre.o: mini_fhtre.c iio.h
minimize.o: minimize.c fail.c minimize_gsl.c smapa.h
mnehs.o: mnehs.c iio.h xmalloc.c fail.c getpixel.c bicubic.c rpc.c \
 xfopen.c smapa.h parsenumbers.c pickopt.c
mnehs_ms.o: mnehs_ms.c iio.h xmalloc.c fail.c getpixel.c bicubic.c \
 smapa.h parsenumbers.c pickopt.c
morsi_demo.o: morsi_demo.c iio.h
morsi_demo2.o: morsi_demo2.c iio.h
ofc.o: ofc.c iio.h fragments.c vvector.h
overflow.o: overflow.c iio.h fail.c xmalloc.c xfopen.c
overpoints.o: overpoints.c iio.h fail.c xfopen.c xmalloc.c smapa.h
pairhom.o: pairhom.c iio.h xmalloc.c fail.c parsenumbers.c pickopt.c
pairsinp.o: pairsinp.c iio.h xmalloc.c fail.c parsenumbers.c
pamle_rec.o: pamle_rec.c smapa.h iio.h
paraflow.o: paraflow.c iio.h fail.c statistics.c xmalloc.c synflow_core.c \
 getpixel.c marching_interpolation.c vvector.h homographies.c smapa.h
percomp.o: percomp.c fft.c iio.h fail.c xmalloc.c
periodize.o: periodize.c iio.h xmalloc.c fail.c
perms.o: perms.c
pickopt.o: pickopt.c
plyflatten.o: plyflatten.c xmalloc.c fail.c smapa.h iio.h
plyroads.o: plyroads.c fail.c xfopen.c parsenumbers.c xmalloc.c
plyroads_mini.o: plyroads_mini.c parsenumbers.c xmalloc.c fail.c
pmba.o: pmba.c xfopen.c fail.c parsenumbers.c xmalloc.c pickopt.c
pmba2.o: pmba2.c xfopen.c fail.c parsenumbers.c xmalloc.c pickopt.c
poisson_rec.o: poisson_rec.c iio.h
poisson_recsep.o: poisson_recsep.c iio.h
polygonify.o: polygonify.c ccproc.c abstract_dsf.c xmalloc.c fail.c \
 smapa.h iio.h
posmax.o: posmax.c iio.h
ppsmooth.o: ppsmooth.c iio.h pickopt.c
qlambda.o: qlambda.c iio.h fragments.c getpixel.c
r.o: r.c iio.h
raddots.o: raddots.c iio.h
radphar.o: radphar.c iio.h xmalloc.c fail.c
rancloud.o: rancloud.c iio.h fragments.c
ranrecs.o: ranrecs.c iio.h fragments.c
really_simplest_inpainting.o: really_simplest_inpainting.c iio.h
remove_small_cc.o: remove_small_cc.c abstract_dsf.c iio.h pickopt.c
replicate.o: replicate.c iio.h xmalloc.c fail.c
rgfield.o: rgfield.c iio.h xmalloc.c fail.c smapa.h
rgfields.o: rgfields.c iio.h xmalloc.c fail.c smapa.h
rgfieldst.o: rgfieldst.c iio.h xmalloc.c fail.c smapa.h
rpc.o: rpc.c xfopen.c fail.c smapa.h
rpc_angpair.o: rpc_angpair.c rpc.c xfopen.c fail.c smapa.h xmalloc.c \
 iio.h
rpc_curvpair.o: rpc_curvpair.c rpc.c xfopen.c fail.c smapa.h xmalloc.c \
 iio.h
rpc_epicyl.o: rpc_epicyl.c rpc.c xfopen.c fail.c smapa.h
rpc_errfilt.o: rpc_errfilt.c parsenumbers.c xmalloc.c fail.c rpc.c \
 xfopen.c smapa.h pickopt.c
rpc_errpair.o: rpc_errpair.c rpc.c xfopen.c fail.c smapa.h xmalloc.c \
 iio.h
rpc_errsingle.o: rpc_errsingle.c rpc.c xfopen.c fail.c smapa.h xmalloc.c \
 iio.h
rpc_eval.o: rpc_eval.c rpc.c xfopen.c fail.c smapa.h
rpc_mnehs.o: rpc_mnehs.c iio.h getpixel.c bicubic.c rpc.c xfopen.c fail.c \
 smapa.h tiffu.c xmalloc.c
rpc_pm.o: rpc_pm.c iio.h getpixel.c bicubic.c rpc.c xfopen.c fail.c \
 smapa.h tiffu.c xmalloc.c
rpc_pmn.o: rpc_pmn.c rpc.c xfopen.c fail.c smapa.h tiffu.c xmalloc.c \
 iio.h
rpc_warpab.o: rpc_warpab.c getpixel.c bicubic.c rpc.c xfopen.c fail.c \
 smapa.h iio.h xmalloc.c
rpc_warpabt.o: rpc_warpabt.c getpixel.c bicubic.c rpc.c xfopen.c fail.c \
 smapa.h tiffu.c iio.h xmalloc.c
rpcflow.o: rpcflow.c iio.h xmalloc.c fail.c rpc.c xfopen.c smapa.h
rpchfilt.o: rpchfilt.c parsenumbers.c xmalloc.c fail.c rpc.c xfopen.c \
 smapa.h
rpcparcheck.o: rpcparcheck.c xmalloc.c fail.c parsenumbers.c rpc.c \
 xfopen.c smapa.h
rpctest.o: rpctest.c iio.h xmalloc.c fail.c rpc.c xfopen.c smapa.h
satproj.o: satproj.c fail.c getpixel.c iio.h xmalloc.c parsenumbers.c \
 pickopt.c
scheme_plap.o: scheme_plap.c fail.c iio.h xmalloc.c
setdim.o: setdim.c iio.h
sfblur.o: sfblur.c xmalloc.c fail.c iio.h
shadowcast.o: shadowcast.c iio.h pickopt.c
sheartilter.o: sheartilter.c
shuntingyard.o: shuntingyard.c fail.c xmalloc.c
simplest_inpainting.o: simplest_inpainting.c iio.h
sort.o: sort.c iio.h
sphereheights.o: sphereheights.c iio.h fail.c xmalloc.c
spline_double.o: spline_double.c
starfield.o: starfield.c random.c iio.h xmalloc.c fail.c
startracker.o: startracker.c xmalloc.c fail.c getpixel.c ok_list.c grid.c \
 iio.h
sterint.o: sterint.c iio.h xmalloc.c fail.c
strt.o: strt.c xmalloc.c fail.c iio.h pickopt.c
sum.o: sum.c iio.h
tcregistration.o: tcregistration.c iio.h
tdip.o: tdip.c iio.h strt.c xmalloc.c fail.c smapa.h random.c pickopt.c
testgblur.o: testgblur.c gblur.c iio.h fail.c xmalloc.c vvector.h
thinpa_recsep.o: thinpa_recsep.c smapa.h iio.h
tiff_octaves_notest.o: tiff_octaves_notest.c
tiff_octaves_old.o: tiff_octaves_old.c
tiffu.o: tiffu.c
tiloct.o: tiloct.c
tilt_and_shear.o: tilt_and_shear.c bicubic.c getpixel.c fail.c xmalloc.c \
 iio.h
tiny_lure.o: tiny_lure.c getpixel.c blur.c fail.c xmalloc.c smapa.h iio.h
tregistration.o: tregistration.c iio.h
tvint.o: tvint.c smapa.h iio.h pickopt.c
unalpha.o: unalpha.c iio.h
uncrop.o: uncrop.c iio.h
unshadow.o: unshadow.c iio.h pickopt.c
utm.o: utm.c
vecov_lm.o: vecov_lm.c iio.h fail.c xmalloc.c pickopt.c
vecstack.o: vecstack.c iio.h
watermark.o: watermark.c xmalloc.c fail.c xfopen.c font_6x12.c iio.h
wgs84.o: wgs84.c
zeropad.o: zeropad.c xmalloc.c fail.c iio.h
zoombil.o: zoombil.c fail.c xmalloc.c bilinear_interpolation.c smapa.h \
 iio.h
