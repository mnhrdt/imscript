__shadowcast = 0  # globally accessible C function

def __setup_functions():
	global __shadowcast
	if __shadowcast != 0: return

	from os.path import dirname
	from ctypes import CDLL, c_float, c_int, c_char_p, c_void_p
	from numpy.ctypeslib import ndpointer
	from sys import platform

	e = "so" if platform != "darwin" else "dylib"
	L = CDLL(f"{dirname(__file__)}/cshadowcast.{e}")
	S = L.cast_shadows
	S.argtypes = [ndpointer(c_float), c_int,c_int, c_float,c_float,c_float]
	S.restype = None
	__shadowcast = S

def shadowcast(
		x,  # input DEM (numpy array)
		α,  # sun azimuth, in degrees (0=north=up, 90=east=right, etc)
		θ   # sun elevation, in degrees (0=horizon, 90=zenith)
		):

	__setup_functions()
	from numpy import ascontiguousarray
	from numpy.ctypeslib import as_array

	p = 0
	q = 0
	r = 1
	if θ >= 0 and θ < 90:
		from math import sin, cos, radians
		p = cos(radians(θ)) * cos(radians(α))
		q = cos(radians(θ)) * sin(radians(α))
		r = sin(radians(θ))

	assert len(x.shape) == 2

	h = x.shape[0]
	w = x.shape[1]
	X = ascontiguousarray(x, dtype='float32')
	__shadowcast(X, w, h, p, q, r)
	y = as_array(X, (h,w)).copy()
	return y

def __export_shadowcast():
	import sys
	sys.modules[__name__] = shadowcast

version = 1
__export_shadowcast()
