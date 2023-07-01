
def programs():
	from subprocess import check_output as sh
	return sh("make -np|grep BIN\ :|cut -c8-", shell=True, text=True)


from setuptools.command.install import install
class myinstall(install):
	def run(self):
		import os
		os.system("mkdir -p bin")
		os.system(f"make -j {programs()}")
		install.run(self)


import setuptools
setuptools.setup(
	name = "imscript",
	version = "0.8",
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
	data_files = [('bin', programs().split())]
)
