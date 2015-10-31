#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include "xfopen.c"
#include "parsenumbers.c"

static int getlinen(char *l, int n, FILE *f)
{
	int c, i = 0;
	while (i < n-1 && (c = fgetc(f)) != EOF && c != '\n')
		l[i++] = c;
	l[i] = '\0';
	return i;
}

static void apply_transform(double y[3], double A[12], double x[3])
{
	y[0] = (A[0]*x[0] + A[1]*x[1] + A[2]*x[2] )+ A[3];
	y[1] = (A[4]*x[0] + A[5]*x[1] + A[6]*x[2] )+ A[7];
	y[2] = (A[8]*x[0] + A[9]*x[1] + A[10]*x[2])+ A[11];
}

struct ply_property {
	enum {UCHAR,FLOAT,DOUBLE,UNKNOWN} type;
	char name[0x100];
	size_t len;
};

static bool parse_property_line(struct ply_property *t, char *buf)
{
	char typename[0x100];
	bool r = 2 == sscanf(buf, "property %s %s\n", typename, t->name);
	t->type = UNKNOWN; t->len = 0;
	if (0 == strcmp(typename, "uchar")) { t->type = UCHAR;  t->len = 1;}
	if (0 == strcmp(typename, "float")) { t->type = FLOAT;  t->len = 4;}
	if (0 == strcmp(typename, "double")){ t->type = DOUBLE; t->len = 8;}
	if (0 == strcmp(typename, "list")){ r = 0; } // skip the face properties
	return r;
}


// reads buf and updates the ply header parameters and properties
// 	 the utm zone
// 	 the vertex properties stored in *t (numnber of properties nprop), 
// 	   *t contains the names and sizes of each field (in bytes)
// 	 the number of vertexes
// 	 isbin if the encoding is binary
// 	 is_ply is set only after the end_header is found
static void process_a_ply_header_line(char *buf, char* utm, int *nvertex, int *isbin, int *is_ply, int *nprop, struct ply_property *t)
{
	if (0 == strcmp(buf, "format binary_little_endian 1.0")) *isbin=1;
	else if (0 == strcmp(buf, "format ascii 1.0")) *isbin=0;
	else if (parse_property_line(t+*nprop, buf)) *nprop += 1;
	else if (0 == strncmp(buf, "comment projection:", 19)) 
		sscanf(buf, "comment projection: UTM %s", utm);
	else if (buf==strstr(buf, "element vertex "))
		sscanf(buf, "element vertex %d", nvertex);

	if (0 == strcmp(buf, "end_header")) *is_ply = 1;
}




int get_record(FILE *f_in, int isbin, struct ply_property *t, int n, double *data)
{
	int rec = 0;
	if(isbin) {
		for (int i = 0; i < n; i++) {
			switch(t[i].type) {
				case UCHAR: {
						    unsigned char X;
						    rec += fread(&X, 1, 1, f_in);
						    data[i] = X;
						    break; }
				case FLOAT: {
						    float X;
						    rec += fread(&X, sizeof(float), 1, f_in);
						    data[i] = X;
						    break; }
				case DOUBLE: {
						     double X;
						     rec += fread(&X, sizeof(double), 1, f_in);
						     data[i] = X;
						     break; }
				default: { break; }
			}
		}
	} else {
		int i=0;
		while (i < n && !feof(f_in)) {
			rec += fscanf(f_in,"%lf", &data[i]);  i++;
		}
	}
	return rec;
}




void put_record(FILE *f_out, int isbin, struct ply_property *t, int n, double *data)
{
	for (int i = 0; i < n; i++) {
		switch(t[i].type) {
			case UCHAR: {
					    unsigned char X = data[i];
					    if (isbin) fwrite(&X, 1, 1, f_out);
					    else fprintf(f_out,"%d", X);
					    break; }
			case FLOAT: {
					    float X = data[i];
					    if (isbin) fwrite(&X, sizeof(float), 1, f_out);
					    else fprintf(f_out,"%a", X);
					    break; }
			case DOUBLE: {
					     double X = data[i];
					     if (isbin) fwrite(&X, sizeof(double), 1, f_out);
					     else fprintf(f_out,"%a", X);
					     break; }
			default: { break; }
		}
		if (!isbin && i<n-1) fprintf(f_out," ");
	}
	if (!isbin) fprintf(f_out,"\n");
}




int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s \"a1 .. a12\" [in [out]]\n", *v);
		//                         0   1            2   3
		return 1;
	}
	char *affstring    = v[1];
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	double A[12];
	read_n_doubles_from_string(A, affstring, 12);
	fprintf(stderr, "3d motion:\n\t%g %g %g %g\n"
			"\t%g %g %g %g\n\t%g %g %g %g\n",
			A[0], A[1], A[2], A[3], A[4], A[5],
			A[6], A[7], A[8], A[9], A[10], A[11]);

	FILE *fi = xfopen(filename_in, "r");
	FILE *fo = xfopen(filename_out, "w");

	int n, nmax = 10000;
	char buf[nmax];
	int cx = 0, maxconvert = 0;


	// try to determine if we are dealing with a ply
	char utm[3];
	int is_ply = 0;
	int isbin = 0;
	int nprop = 0;
	struct ply_property t[100];


	// read the ply header or the whole file if not binary encoded or not ply at all
	while ((n = getlinen(buf, nmax, fi)))
	{

		process_a_ply_header_line(buf, utm, &maxconvert, &isbin, &is_ply, &nprop, t);
		//    the following lines process a register if formatted adequately
		//    (x,y,z, ...) 
		//    otherwise print the buffer (the line) in the output
		double x[3], y[3];
		int d, r = sscanf(buf, "%lf %lf %lf %n", x, x+1, x+2, &d);
		if (r == 3 && (!maxconvert || cx < maxconvert)) {
			apply_transform(y, A, x);
			fprintf(fo, "%lf %lf %lf ", y[0], y[1], y[2]);
			fputs(buf + d, fo);
			cx += 1;
		} else
			fputs(buf, fo);
		fputc('\n', fo);

		// is_ply is set to true at end_header
		// if it's a binary ply we need different code
		if(is_ply && isbin) break;
	}


	// continue processing the binary part of the ply
	if(is_ply && isbin) {
		double data[nprop];

		while ( cx < maxconvert && 
				nprop == get_record(fi, isbin, t, nprop, data) ) {
			double y[3];
			apply_transform(y, A, data);
			data[0]=y[0];
			data[1]=y[1];
			data[2]=y[2];
			put_record(fo, isbin, t, nprop, data);
			cx++;
		}

		// continue copying the faces of the ply (if present)
		unsigned char X;
		while (fread (&X, 1, 1, fi)) 
			fwrite(&X, 1, 1, fo);
	}

	fprintf(stderr, "converted %d points\n", cx);
	if (maxconvert)
		assert(maxconvert == cx);


	xfclose(fo);
	xfclose(fi);
	return 0;
}
