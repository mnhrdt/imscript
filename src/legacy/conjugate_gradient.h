// the type of a function f that evaluates a linear map y=f(x)
typedef void (*linear_map_t)(double *y, double *x, int n, void *e);

// find x such that Ax=b
void conjugate_gradient(double *x, linear_map_t A, double *b, int n, void *e);

// fancier version with initialization and control parameters
void fancy_conjugate_gradient(double *x,
		linear_map_t A, double *b, int n, void *e,
		double *x0, int max_iter, double min_residual);
