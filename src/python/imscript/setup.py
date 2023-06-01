
from setuptools.command.install import install
class myinstall(install):
	def run(self):
		import os
		os.system("mkdir -p bin")
		os.system("make -j bin/plambda bin/qauto bin/morsi")
		install.run(self)


import setuptools
setuptools.setup(
	name = "imscript",
	version = "0.5",
	author = "Enric Meinhardt-Llopis",
	author_email = "enric.meinhardt@fastmail.com",
	description = "Image Processing with Unix Pipes",
	url = "https://github.com/mnhrdt/imscript",
	classifiers = [
		"Operating System :: OS Independent",
		"License :: OSI Approved :: GNU Affero General Public License v3",
		"Topic :: Scientific/Engineering :: Image Processing",
		"Topic :: Scientific/Engineering :: Mathematics"
		],
	cmdclass = {'install' : myinstall},
	data_files = [('bin', ['bin/plambda', 'bin/qauto', 'bin/morsi'])]
)
