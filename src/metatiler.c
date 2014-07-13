// metatiler: run a iio program tile-wise on a tiled tiff

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tiffio.h>
#include "iio.h"

#define CMDLINE_MAX 10000

#define MARKER_INPUT  '^'
#define MARKER_OUTPUT '@'

static void add_item_to_cmdline(char *cmdline, char *item, char *fileprefix)
{
	//fprintf(stderr, "ADD \"%s%s\"\n", fileprefix?fileprefix:"", item);
	if (*cmdline)
		strncat(cmdline, " ", CMDLINE_MAX);
	if (fileprefix)
		strncat(cmdline, fileprefix, CMDLINE_MAX);
	strncat(cmdline, item, CMDLINE_MAX);
}

static void fill_subs_cmdline(char *cmdline, char *cmd, char *fileprefix,
		char **fns_in, int n_in, char **fns_out, int n_out)
{
	*cmdline = 0;
	char *tok = strtok(cmd, " ");
	do {
		if (*tok=='>' || *tok=='|' || *tok=='<') {
			fprintf(stderr, "ERROR: must be a single "
					"command line\n");
			exit(1);
		} else if (*tok == MARKER_INPUT) {
			int idx = atoi(tok+1)-1;
			if (idx >= 0 && idx < n_in)
			add_item_to_cmdline(cmdline, fns_in[idx], fileprefix);
		} else if (*tok == MARKER_OUTPUT) {
			int idx = atoi(tok+1)-1;
			if (idx >= 0 && idx < n_out)
			add_item_to_cmdline(cmdline, fns_out[idx], fileprefix);
		} else
			add_item_to_cmdline(cmdline, tok, NULL);
	} while (tok = strtok(NULL, " "));
	fprintf(stderr, "CMDLINE = \"%s\"\n", cmdline);
}

static char *create_temporary_directory(void)
{
	return "/tmp/MTD/";
	//return "/tmp/metafilter_temporary_directory/";
}

static void extract_tile(char *tpd, char *filename, int idx)
{
	// 1. open large tiff image "filename"
	// 2. locate idexed tile
	// 3. read indexed tile data
	// 4. close large tiff image
	// 5. create new small tiff image "tpd/filename"
	// 6. write data to new small tiff image
}

static void paste_tile(char *filename, int idx, char *tpd)
{
	// 1. read data from small tiff file
	// 2. open large tiff file
	// 3. locate position of indexed tile
	// 4. rewrite indexed tile with new data
	// 5. save and close the large tiff file
}

void metatiler(char *command, char **filenames_in, int n_in,
		char **filenames_out, int n_out)
{
	// build command line (will be the same at each run)
	char cmdline[CMDLINE_MAX];
	char *tpd = create_temporary_directory();
	fill_subs_cmdline(cmdline, command, tpd,
			filenames_in, n_in, filenames_out, n_out);

	// determine tile geometry
	int ntiles = 4;//=get_number_of_tiles(*filename_in)

	// create output files as empty files with the required tile geometry
	// to solve: what pixel dimension and type?
	// solution: run the program with the first tile and see what happens

	for (int i = 0; i < ntiles; i++)
	{
		// 1. extract ith tile from input images
		// 2. run cmdline
		fprintf(stderr, "run \"%s\"\n", cmdline);
		// 3. paste results to output files
	}
}

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "usage:\n\t"
			"%s \"CMD ^1 ^2 @1\" in1 in2 -- out1\n", *argv);
		//       0   1               2   3   ...
		return 1;
	}

	// get input arguments
	char *command = argv[1];
	int n_in = 0, n_out = 0;
	char *filenames_in[argc];
	char *filenames_out[argc];
	for (int i = 2; i < argc && strcmp(argv[i], "--"); i++)
		filenames_in[n_in++] = argv[i];
	for (int i = 3+n_in; i < argc; i++)
		filenames_out[n_out++] = argv[i];

	// print debug info
	fprintf(stderr, "%d input files:\n", n_in);
	for (int i = 0; i < n_in; i++)
		fprintf(stderr, "\t%s\n", filenames_in[i]);
	fprintf(stderr, "%d output files:\n", n_out);
	for (int i = 0; i < n_out; i++)
		fprintf(stderr, "\t%s\n", filenames_out[i]);
	fprintf(stderr, "COMMAND = \"%s\"\n", command);

	// run program
	metatiler(command, filenames_in, n_in, filenames_out, n_out);

	// exit
	return 0;
}
