#include <math.h>

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
		int datadim,       // dimension of each data point
		int n              // number of data points
		int modeldim,      // number of model parameters
		model_evaluating_function *mev,
		int nfit,          // data points needed to produce the model

		// decoration
		void *user
		);


static int random_index(int a, int b)
{
	int r = a + rand()%(b - a + 1);
	assert(r >= a);
	assert(r < b);
	return r;
}

static bool are_different(int *t, int n)
{
	qsort(t, n, sizeof*t, compare_ints);
}

static void fill_random_indices(int *idx, int ni, int a, int b)
{
	int safecount = 0;
	while (safecount < 10 && !are_different(idx, ni))
	{
		for (int i = 0; i < ni; i++)
			idx[i] = random_index(a, b);
		safecount += 1;
	}
	if (safecount == 10)
		fail("could not generate any model");
}


// RANSAC
//
// Given a list of data points, find the parameters of a model that fits to
// those points.  Several models are tried, and the model with the highest
// number of inliers is kept.
bool ransac(
		// output
		float *out_error,  // error measure of the set of inliers
		int *out_ninliers, // number of inliers
		bool *out_mask,    // array mask identifying the inliers

		// input data
		float *data,       // array of input data

		// input context
		int datadim,       // dimension of each data point
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
		)
{
	int best_ninliers = 0;
	float best_model[modeldim];
	bool *best_mask = xmalloc(n * sizeof*best_mask);

	for (int i = 0; i < ntrials; i++)
	{
		int indices[nfit];
		// TODO fisher yates shuffle and traverse it by blocks of nfit
		fill_random_indices(indices, nfit, 0, n);

		float x[nfit*datadim];
		for (int j = 0; j < nfit; j++)
		for (int k = 0; k < datadim; k++)
			x[datadim*j + k] = data[datadim*indices[j] + k];

		float model[modeldim];
		mgen(model, modeldim, x, datadim, nfit, user);

		float this_error;
		int n_inliers;
		for (int j = 0; j < n; j++)
			out_mask[j] = false;
		ransac_trial(&this_error, &n_inliers, best_mask,
				data, model, datadim, n, modeldim, mev, nfit,
				user);

		if (n_inliers > best_inliers)
		{
		}


	}
}





