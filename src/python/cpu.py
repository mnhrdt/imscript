# cpu(x) displays the image x using the cpu image viewer
#
# TODO 0. add an option to make cpu call non-blocking
# TODO 1. make it work for several images: cpu([a,b,c]) flips them
# TODO 2. traverse the call stack to recover the variable names of each image
#         use inspect.stack() instead of get_ipython()
# TODO 3. v = cpu(x) returns the crop and the contrast change
#         v is a "view", determining a crop in space and color, for example
#         v = (cmin, cmax, xmin, xmax, ymin, ymax)
# TODO 4. vendor the static executable (probably αpε) of cpu


# internal function to do some notebook magic (only for gray images, by now)
def __heuristic_reshape(s):
	try:
		w = get_ipython().all_ns_refs[0]['w']
		h = get_ipython().all_ns_refs[0]['h']
		if s[0] == w*h:
			return (h,w)
		else:
			return s
	except (NameError, KeyError):
		return s

# API
def cpu(x):
	if x == "version":
		global version
		return version

	if x.shape[0] > 4000:
		x = x.reshape(__heuristic_reshape(x.shape))

	import tempfile, iio, os

	f = f"{tempfile.NamedTemporaryFile().name}.tiff"
	c = f"cpu {f} 2>/dev/null ; rm {f}"

	iio.write(f, x)
	os.system(c)

# auxiliary function to import the module as a function
def __export_cpu():
	import sys
	sys.modules[__name__] = cpu

version = 1
__export_cpu()
