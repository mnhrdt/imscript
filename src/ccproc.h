int ccproc(
		int *out_size,          // total size of each CC
		int *out_bdsize,        // boundary size of each CC
		int *out_all,           // array of all the indexes, connected
		int *out_first,         // array of first indexes of each CC
		int *out_idx,           // image with the indices of each region
		float *x, int w, int h, // input image
		int (eq)(float,float)   // input equivalence relation
	);                              // return value = number of components

int bfollow(
		int *out_boundary,      // (x,y) coordinates along boundary
		float *x, int w, int h, // input image
		int (eq)(float,float),  // input equivalence relation
		int i, int j            // pixel inside the CC to follow
	);                              // return value = length of boundary
