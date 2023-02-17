def cpu(x):
	import tempfile, iio, os

	f = f"{tempfile.NamedTemporaryFile().name}.tiff"
	c = f"cpu {f} 2>/dev/null ; rm {f}"

	iio.write(f, x)
	os.system(c)

import sys
sys.modules[__name__] = cpu
