// GENERAL NOTATION
// int w, h : dimensions of the image domain
// float *g : array of length 3*w*h containin the w*h metric tensors

// scalar product of two vectors at a point
//	g: image of size 3*w*h containing the metric tensor
//	w: image width
//	h: image heigth
//	x: base point
//	p: first vector
//	q: second vector
float rip_product(float *g, int w, int h, float x[2], float p[2], float q[2]);

// length of a sampled curve
//	g: image of size 3*w*h containing the metric tensor
//	w: image width
//	h: image heigth
//	x: array of coordinates of the successive points in the curve
//	n: number of points
float rip_curve_length(float *g, int w, int h, float *x, int n);

// energy of a sampled curve
//	g: image of size 3*w*h containing the metric tensor
//	w: image width
//	h: image heigth
//	x: array of coordinates of the successive points in the curve
//	n: number of points
float rip_curve_energy(float *g, int w, int h, float *x, int n);

// next point in a polygonal geodesic (follow Snell's law)
//	g: image of size 3*w*h containing the metric tensor
//	w: image width
//	h: image heigth
//	c: index of the cell at point x
//	x: base point
//	v: direction from x
//	out_x: output point
//	out_v: direction from the output point
//	return value: index of the cell at point "out_x"
int rip_next_snell(float *g, int w, int h, int c, float x[2], float v[2],
		float out_x[2], float out_v[2]);

// exponential map
//	g: image of size 3*w*h containing the metric tensor
//	w: image width
//	h: image heigth
//	c: index of the cell at point x
//	x: base point
//	v: speed from x
int rip_exponential(float *g, int w, int h, int c, float x[2], float v[2],
		float 
		);

float rip_logarithm(float *g, int w, int h, int c, float p[2], float q[2],
		float out_v[2]);
