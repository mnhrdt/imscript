def cpu(x):
	import tempfile, piio, os

	f = f"{tempfile.NamedTemporaryFile().name}.tiff"
	c = f"cpu {f} 2>/dev/null ; rm {f}"

	piio.write(f, x)
	os.system(c)
