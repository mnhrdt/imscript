// generic function
// evaluate the error of a datapoint according to a model
void ransac_error_evaluation(
		float *model,
		int modeldim,
		float *datapoint,
		int datadim,
		void *user
		);


// generic function
// evaluate the error of a model
void ransac_model_evaluating_function(
		float a,
		);

// generic function
// compute the model defined from a few data points
void ransac_model_generating_function(
		float *out_model,  // parameters of the computed model
		int modeldim,      // number of parameters of the model (unused)
		float *data,       // data points
		int datadim,       // dimension of data points (unused)
		int nfit,          // number of data points (probably unused)
		void *user
		);


void ransac_trial(
		// output
		float *out_error,  // error measure of the set of inliers
		int *out_ninliers, // number of inliers
		bool *out_mask,    // array mask identifying the inliers

		// input data
		float *data,       // array of input data
		float *model,      // parameters of the model

		// input context
		int dim,           // dimension of each data point
		int n              // number of data points
		int modeldim,      // number of model parameters
		model_evaluating_function *mev,
		int nfit,          // data points needed to produce the model

		// decoration
		void *user
		);



// RANSAC
//
// Given a list of data points, find the parameters of a model that fits to
// those points.  Several models are tried, and the model with least error is
// kept.
bool ransac(
		// output
		float *out_error,  // error measure of the set of inliers
		int *out_ninliers, // number of inliers
		bool *out_mask,    // array mask identifying the inliers

		// input data
		float *data,       // array of input data

		// input context
		int dim,           // dimension of each data point
		int n,             // number of data points
		int modeldim,      // number of model parameters
		ransac_model_evaluating_function *mev,
		ransac_model_generating_function *mgen,
		int nfit,          // data points needed to produce a model

		// input parameters
		int ntrials,       // number of models to try
		int min_inliers,   // minimum allowed number of inliers
		float max_error,   // maximum allowed error

		// decoration
		ransac_model_accepting_function *macc,
		void *user
		);





