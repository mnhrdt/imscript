
from setuptools.command.build_ext import build_ext
class my_build(build_ext):
	def run(self):
		N = "shadowcast"                      # name
		O = "-O3 -march=native -shared -fPIC" # compilation options
		E = "so"                              # .ext of shared objects
		import sys, os
		if sys.platform == "darwin":
			E = "dylib"
			O = f"{O} -undefined dynamic_lookup"
		os.system(f"cc {O} {N}.c -o c{N}.{E}")
		os.system(f"cp c{N}.{E} {self.build_lib}/c{N}.{E}")

from setuptools import Extension, setup
setup(
	name="shadowcast",
	version='1',
	author="Enric Meinhardt-Llopis",
	author_email="enric.meinhardt@ens-paris-saclay.fr",
	description="python wrapper to shadowcast",
	url='https://github.com/mnhrdt/imscript',
	classifiers=[
		"Operating System :: OS Independent",
	],
	py_modules=["shadowcast"],
	ext_modules=[ Extension(name="shadowcast", sources=["shadowcast.c"]) ],
	cmdclass={'build_ext': my_build}
)

