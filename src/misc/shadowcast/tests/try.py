import shadowcast, iio, cpu, os
n = "fuji.npy"
if not os.path.exists(n):
	iio.write(n, iio.read("http://gabarro.org/img/fuji.tif"))
x = iio.read(n)[:,:,0]
f = 90/4000
y = shadowcast(x*f, 45, 45)/f
cpu(y)
