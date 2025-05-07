def extensions():
	from setuptools import Extension
	return [Extension("libshadowcast", ["shadowcast.c"])]


from distutils.command.install_lib import install_lib as _install_lib

def batch_rename(src, dst, src_dir_fd=None, dst_dir_fd=None):
	'''Same as os.rename, but returns the renaming result.'''
	import os
	os.rename(src, dst,
		src_dir_fd=src_dir_fd,
		dst_dir_fd=dst_dir_fd)
	return dst

class _CommandInstall(_install_lib):
	def __init__(self, *args, **kwargs):
		_install_lib.__init__(self, *args, **kwargs)

	def install(self):
		# let the distutils' install_lib do the hard work
		outfiles = _install_lib.install(self)

		# batch rename the outfiles:
		# for each file, match string between
		# second last and last dot and trim it
		import re
		matcher = re.compile('\.([^.]+)\.so$')
		return [batch_rename(file, re.sub(matcher, '.so', file))
			for file in outfiles]


from setuptools import setup
setup(name="shadowcast",
	version='1',
	author="Enric Meinhardt-Llopis",
	author_email="enric.meinhardt@fastmail.com",
	description="Python wrapper to shadowcast",
	url='https://github.com/mnhrdt/imscript',
	classifiers=[
		"Operating System :: OS Independent",
	],
	py_modules=['shadowcast'],
	ext_modules=extensions(),
	cmdclass={
		'install_lib': _CommandInstall,
	},
)

