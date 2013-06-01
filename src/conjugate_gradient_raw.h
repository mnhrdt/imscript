// the type of a function f that evaluates a linear map y=f(x)
typedef void (*linear_map)(double *y, double *x, int n);

// find x such that Ax=b
void conjugate_gradient(double *x, linear_map A, double *b, int n);
