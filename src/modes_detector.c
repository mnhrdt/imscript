#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "smapa.h"
SMART_PARAMETER(BINWIDTH, 5)
SMART_PARAMETER(VERBOSE, 0)

// compute the cumulative sum of an array of ints
static void cumsum(int *output, int *input, int n)
{
    output[0] = input[0];
    for (int i = 1; i < n; i++)
        output[i] = input[i] + output[i-1];
}


// compute the Kullback-Leibler divergence (also called relative entropy)
static float entropy(float r, float p)
{
    if (p == 0) return 0;
    if (p == 1) return 0;
    if (r == 0) return -log(1-p);
    if (r == 1) return -log(p);
    return r * log(r/p) + (1-r) * log((1-r)/(1-p));
}


// compute the binomial tail B(M, k, p) TODO
static long double binomial_tail(int M, int k)
{
    long double s = 0;
    for (int i = k; i<=M; i++) {
        s += 0;
    }
    return s;
}


// Function running over all intervals contained in [1, n],
// computing if they are meaningful intervals or gaps.
// The results are stored in the matrices intervals and entropy. In the matrix
// intervals, we use the following convention :
// -meaningful interval : 2
// -meaningful gap : 1
// -neither meaningful interval or gap : 0
static void browse_intervals(int *histo, int n, float epsilon, char **intervals, float **entropies)
{
    // cumulative histogram
    int *cumhist = malloc(n*sizeof(*cumhist));
    cumsum(cumhist, histo, n);

    // total number of samples in the histogram
    int nsamples = cumhist[n - 1];

    float thresh = logf((float) n * (n+1) / (2 * epsilon)) / nsamples;
    if (VERBOSE()) fprintf(stderr, "thresh: %f\n", thresh);

    for (int a = 0; a < n; a++)
    for (int b = a; b < n; b++) {
        int k = cumhist[b] - cumhist[a] + histo[a];
        float r = (float) k / nsamples;
        float p = (float) (1 + b - a) / n;
        float e = entropy(r, p);
        //if (VERBOSE()) fprintf(stderr, "k, r, p, e: %d, %f %f %f\n", k, r, p, e);

        // save the entropy value e
        entropies[a][b] = e;

        // memorize if the interval [a, b] is a meaningful interval or gap
        if (e >= thresh) {
            if (r >= p) intervals[a][b] = 2;
            else intervals[a][b] = 1;
        }
        else
            intervals[a][b] = 0;
    }
}


// This function "propagates" the 1 values such that if an interval contains a gap, its
// marker in the "intervals" matrix goes to 0
static void spread_gaps(int n, char **intervals)
{
    for (int a = 0; a < n; a++)
    for (int b = a; b < n; b++) {
        if (intervals[a][b] == 1)  // it's a gap
            for (int i = a; i >= 0; i--)
            for (int j = b; j < n; j++)
                intervals[i][j] = 0;
    }
}


// For each mode (marker = 2), this function looks at the intervals inside, and
// if there is one with bigger entropy, it discards the
// mode (by putting the marker to 0)
static void discard_modes(int n, char **intervals, float **entropies)
{
    for (int a = 0; a < n; a++)
    for (int b = a; b < n; b++)
        if (intervals[a][b] == 2) {  // [a, b] is a possible maximal mode
            float e = entropies[a][b];
            for (int i = a; i <= b; i++)
            for (int j = b; j >= i; j--) {
                // if [i, j] is a sub-mode strictly included in [a,b]
                if ((j-i < b-a) && (intervals[i][j] == 2)) {
                    if (entropies[i][j] > e) {  // [a,b] is not a maximal mode
                        intervals[a][b] = 0;
                    } else {
                        intervals[i][j] = 0;
                    }
                }
            }
        }
}


struct interval {
    char a;
    char b;
    float nfa;
};


int compare_meaningful_intervals(const void *pa, const void *pb)
{
    const struct interval *a = pa;
    const struct interval *b = pb;
    return a->nfa - b->nfa;
}


// Arguments:
// int **modes: output array of size nmodes_max, 2
// int nmodes_max: maximal number of modes we allow to detect
// int *histogram: input histogram
// int n: length of the input histogram
// float epsilon: parameter of the a contrario model
static int histogram_modes(int modes[][2], float nfas[], int nmodes_max, int *histogram, int n, float epsilon)
{
    // n by n array to remember if interval [a, b] is a meaningful interval or gap
    char **intervals = malloc(n*sizeof(*intervals));
    float **entropies = malloc(n*sizeof(*entropies));
    struct interval detections[n];  // array of detected meaningful intervals
    for (int i = 0; i < n; i++) {
        intervals[i] = malloc(n*sizeof(**intervals));
        entropies[i] = malloc(n*sizeof(**entropies));
    }

    // loop over all intervals to compute entropies
    browse_intervals(histogram, n, epsilon, intervals, entropies);
//    if (VERBOSE()) {
//        fprintf(stderr, "entropies:\n");
//        for (int a = 0; a < n; a++) {
//            for (int b = 0; b < a; b++) fprintf(stderr, "      ");
//            for (int b = a; b < n; b++) fprintf(stderr, "%.03f ", entropies[a][b]);
//            fprintf(stderr, "\n");
//        }
//        fprintf(stderr, "intervals:\n");
//        for (int a = 0; a < n; a++) {
//            for (int b = 0; b < a; b++) fprintf(stderr, "  ");
//            for (int b = a; b < n; b++) fprintf(stderr, "%d ", intervals[a][b]);
//            fprintf(stderr, "\n");
//        }
//    }

    // discard intervals containing gaps
    spread_gaps(n, intervals);
//    if (VERBOSE())
//        for (int a = 0; a < n; a++) {
//            for (int b = 0; b < a; b++) fprintf(stderr, "  ");
//            for (int b = a; b < n; b++) fprintf(stderr, "%d ", intervals[a][b]);
//            fprintf(stderr, "\n");
//        }

    // keep maximal modes
    discard_modes(n, intervals, entropies);
//    if (VERBOSE())
//        for (int a = 0; a < n; a++) {
//            for (int b = 0; b < a; b++) fprintf(stderr, "  ");
//            for (int b = a; b < n; b++) fprintf(stderr, "%d ", intervals[a][b]);
//            fprintf(stderr, "\n");
//        }

    // cumulative histogram
    int *cumhist = malloc(n*sizeof(*cumhist));
    cumsum(cumhist, histogram, n);

    // total number of samples in the histogram
    int nsamples = cumhist[n - 1];

    // list meaningful modes
    int nmodes = 0;
    float *x = malloc(n*sizeof(*x)); // there will never be more than n modes
    for (int a = 0; a < n; a++)
    for (int b = a; b < n; b++)
        if (intervals[a][b] > 1) {
            // approximated NFA
            float nfa = n * (n+1) / 2 * exp(-nsamples * entropies[a][b]);
            detections[nmodes].a = a;
            detections[nmodes].b = b;
            detections[nmodes].nfa = nfa;
            nmodes += 1;

            //int k = histo.sum(a,b);
            //float p = (1+b-a)/((float) L) + (b<a);
            //float nfa = binomial_tail(floor(nsamplesM+0.5),k,p)*histo.get_N();
        }

    // sort detected modes according to nfa
    qsort(detections, nmodes, sizeof(*detections), compare_meaningful_intervals);
    if (nmodes > nmodes_max) nmodes = nmodes_max;
    for (int i = 0; i < nmodes; i++) {
        modes[i][0] = detections[i].a;
        modes[i][1] = detections[i].b;
        nfas[i] = detections[i].nfa;
    }

    for (int i = 0; i < n; i++) {
        free(intervals[i]);
        free(entropies[i]);
    }
    free(intervals);
    free(entropies);
    return nmodes;
}


static void compute_histogram(int *histogram, float bin_width, float *points, int npts)
{
    // determine min and max input points
    float min = +INFINITY;
    float max = -INFINITY;
	for (int i = 0; i < npts; i++) {
        if (points[i] < min)
            min = points[i];
        if (points[i] > max)
            max = points[i];
    }

    // initialize histogram bins to zero
    int nbin = 1 + (max - min) / bin_width;
	for (int i = 0; i < nbin; i++)
       histogram[i] = 0;

    // loop over the points to increment histogram bins
	for (int i = 0; i < npts; i++) {
       int k = (points[i] - min) / bin_width;
       if (k >= nbin) k = nbin - 1;
       histogram[k]++;
    }
}


static float histogram_mode_average(int *histogram, int nbin, int mode[2])
{
    assert(0 <= mode[0]);
    assert(mode[0] <= mode[1]);
    assert(mode[1] < nbin);
    float sum = 0;
    int n = 0;
	for (int i = mode[0]; i <= mode[1]; i++) {
       n += histogram[i];
       sum += i*histogram[i];
    }
    return sum / n;
}


static float average_points_between_bounds(float *x, int n, float min, float max)
{
    float sum = 0;
    int c = 0;
	for (int i = 0; i < n; i++) {
        if ((x[i] >= min) && (x[i] <= max)) {
            sum += x[i];
            c++;
            //if (VERBOSE())
            //    fprintf(stderr, "c, x, sum: %d, %f, %f\n", c, x[i], sum);
        }
    }
    if (c > 0)
        return sum / c;
    else
        return NAN;
}


// OUTPUT float y[ny]: y[0] is the number of detected modes. Next values are the centers of the detected modes
// INPUT float x[nx]: nx points
void acontrario_modes_detector(float *y, int ny, float *x, int nx, float epsilon)
{
    // determine min and max input points to compute the number of bins
    float min = +INFINITY;
    float max = -INFINITY;
	for (int i = 0; i < nx; i++) {
        if (x[i] < min)
            min = x[i];
        if (x[i] > max)
            max = x[i];
    }

    // abort if there are no finite input points
    if (!(isfinite(min) && isfinite(max))) {
        y[0] = 0;
        return;
    }

    int nbin = 1 + (max - min) / BINWIDTH();
    if (VERBOSE()) fprintf(stderr, "min, max, nbin: %.0f %.0f %d\n", min, max, nbin);

    // build the histogram
    int *histogram = malloc(nbin*sizeof(*histogram));
    compute_histogram(histogram, BINWIDTH(), x, nx);

    if (VERBOSE()) {
        fprintf(stderr, "min max %.0f %.0f\n", min, max);
        fprintf(stderr, "histogram: ");
	    for (int i = 0; i < nbin; i++)
            fprintf(stderr, "%d ", histogram[i]);
        fprintf(stderr, "\n");
    }

    // run the histogram modes detector
    int modes[ny-1][2];
    float nfas[ny-1];
    int nmodes = histogram_modes(modes, nfas, ny-1, histogram, nbin, epsilon);
    if (VERBOSE()) {
	    for (int i = 0; i < nmodes; i++)
            fprintf(stderr, "[%d, %d] %f\n", modes[i][0], modes[i][1], nfas[i]);
    }

    // compute the mean of the input points in each mode
    y[0] = nmodes;
    for (int i = 0; i < nmodes; i++) {
        //float x = histogram_mode_average(histogram, nbin, modes[i]);
        //y[i+1] = min + (x+.5) * (max - min) / nbin;
        //if (VERBOSE()) fprintf(stderr, "%.2f %.2f\n", x, y[i+1]);
        float xa = min + BINWIDTH() * modes[i][0];
        float xb = min + BINWIDTH() * (modes[i][1] + 1);
        y[i+1] = average_points_between_bounds(x, nx, xa, xb);
        if (VERBOSE()) fprintf(stderr, "mode %.2f, nfa %.3f\n", y[i+1], nfas[i]);
    }
}


#ifdef TEST_MAIN
#include "pickopt.c"

static void help(char *v)
{
    fprintf(stderr, "usage: BINWIDTH=1 %s [-e epsilon] [-m max_nb_modes] v1 ...\n", v);
}


int main(int c, char *v[])
{
    if (c < 2) {
        help(*v);
        return EXIT_FAILURE;
    }

    // parse cmdline arguments
	float epsilon = atof(pick_option(&c, &v, "e", "1"));
	int nmodes_max = atoi(pick_option(&c, &v, "m", "20"));
	//float bin_width = pick_option(&c, &v, "w", 1);

    // read input values
	int npts = c - 1;
	float points[npts];
	for (int i = 0; i < npts; i++)
		points[i] = atof(v[i+1]);

    if (VERBOSE()) {
        fprintf(stderr, "input points: ");
	    for (int i = 0; i < npts; i++)
            fprintf(stderr, "%.0f ", points[i]);
        fprintf(stderr, "\n");
    }

    float *modes = malloc((1+nmodes_max)*sizeof(*modes));
    acontrario_modes_detector(modes, 1+nmodes_max, points, npts, epsilon);
    fprintf(stderr, "found %d modes\n", (int) modes[0]);
    for (int i = 0; i < modes[0]; i++)
        fprintf(stderr, "\t %.2f ", modes[i+1]);
    fprintf(stderr, "\n");
}
#endif
