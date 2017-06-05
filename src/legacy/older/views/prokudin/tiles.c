#include <stdio.h>
#include <stdlib.h>

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s w h nw nh "
		//                          0 1 2 3  4
				"| while read cx x0 y0 xf yf ; do ...\n", *v);
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	int nw = atoi(v[3]);
	int nh = atoi(v[4]);

	int sizex = 2 * w / (nw+1);
	int sizey = 2 * h / (nh+1);
	int offx = sizex / 2;
	int offy = sizey / 2;

	int cx = 0;
	for (int j = 0; j < nh; j++)
	for (int i = 0; i < nw; i++)
	{
		int from[2] = { i * offx, j * offy };
		int to[2] = { from[0] + sizex, from[1] + sizey };
		printf("%d %d %d %d %d\n", cx, from[0], from[1], to[0], to[1]);
		cx += 1;
	}
	return 0;
}
