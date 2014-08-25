#include "fancy_image.h"

int main(void)
{
	char *filename = "/tmp/zcaca.tiff";
	double factor = 2;

       	struct fancy_image *f = fancy_image_open(filename, "rw");

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < f->pd; l++)
	{
		double x = fancy_image_getsample(f, i, j, l);
		x = x * factor;
		fancy_image_setsample(f, i, j, l, x);
	}

	fancy_image_close(f);

	return 0;
}

//// option 1 (return the struct as a local variable)
//struct caca x = caca_open("file");
//caca_process(&x);
//caca_close(&x);
//
//// option 2 (return a pointer to an allocated struct, to be freed later)
//struct caca *x = caca_open("file");
//caca_process(x);
//caca_close(x);
//
//// option 3 (fill-in the fields of the given struct)
//struct caca x;
//caca_open(&x, "file");
//caca_process(&x);
//caca_close(&x);
//
//// option 4 (like 3, but different variable names)
//struct caca x[1];
//caca_open(x, "file");
//caca_process(x);
//caca_close(x);
//
//// option 5 (like 1, but different variable names)
//struct caca x[1] = {caca_open("file")};
//caca_process(x);
//caca_close(x);
//
//// option 6 (like 1, but different variable names)
//struct caca x[1];
//x[0] = caca_open("file");
//caca_process(x);
//caca_close(x);
