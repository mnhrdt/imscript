#include <stdbool.h>


struct lk_parameters {
	int window_radius;     // integer radius of the patch around each point
	float window_variance; // variance of the gaussian weights
	int model_degree;      // degree of local polynomial model
	bool stereo;           // force horizontal displacements

	bool rafa;
	float rafa_sigma;

	bool interpolate_well;
	// whatever else
};


// INPUT DATA:
//
//          in_a: first gray-level input image
//          in_b: second gray-level input image
//             w: width of each input image
//             h: height of each input image
//     in_models: optional model initialization
//    in_nmodels: dimension of model initialization data
//   out_nmodels: desired number of dimensions for model initialization
// lk_parameters: parameters for the Lucas-Kanade algorithm
//
//
// OUTPUT DATA:
//
//  out_models: computed flow models
//
static void lk_core(float *out_models, int out_nmodels
		float *in_a, float *in_b, int w, int h,
		float *in_models, int in_nmodels,
		struct lk_parameters *p)
{
	if (false) { ;
	} else if (p->model_degree = 0 && !p->stereo && !in_models) {
		// optical flow without initialization, locally constant
		//
		// model_0: flow x
		// model_1: flow y
		// model_2: structure tensor trace
		// model_3: structure tensor determinant
		// model_4: structure tensor a
		// model_5: structure tensor b
		// model_6: structure tensor c
		// model_7: right hand side x
		// model_8: right hand side y
		//
		//
		if (out_nmodels == 2)
			lk_classic(out_models, in_a, in_b, w, h,
					p->window_radius, p->window_variance);
	} else if (p->model_degree = 0 && p->stereo && !in_models) {
		// stereo without initialization, locally constant
		//
		// model_0: disparity
		// model_1: structure scalar
		// model_2: right hand side
		//
		//
	} else if (p->model_degree = 1 && !p->stereo && !in_models) {
		// optical flow without initialization, locally affine
		//
		// model_0: flow x
		// model_1: flow y
		// model_2: affinity angle
		// model_3: affinity tilt
		// model_4: affinity lambda
		// model_5: affinity mu
		// model_6: affinity a
		// model_7: affinity b
		// model_8: affinity c
		// model_9: affinity d
		// ...
		//
		//
	} else if (p->model_degree = 1 && p->stereo && !in_models) {
		// stereo without initialization, locally affine
		//
		// etc
	}
}
