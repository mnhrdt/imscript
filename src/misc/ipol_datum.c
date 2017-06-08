
struct ipol_datum {
	enum {
		image, number, string
	} type;
	int w;
	int h;
	int pd;
	float *data;
	float datum;
	char *string;
};

#include <stdlib.h>

void ipol_datum_free(struct ipol_datum *x)
{
	switch(x->type) {
	case image: free(x->data); break;
	case string: free(x->string);
	case number:
		     break;
	}
}
