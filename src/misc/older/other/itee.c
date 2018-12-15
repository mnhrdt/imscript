// itee, copy stdin to stdout, and open an image (using fpan) with this data

#include <stdio.h>
#include <stdlib.h>

int main()
{
	char fname[] = "/tmp/itee_XXXXXX";
	int fd = mkstemp(fname);
	FILE *f = fdopen(fd, "w");
	int c;
	while ((c = getchar()) != EOF)
	{
		putchar(c);
		fputc(c, f);
	}
	putchar(EOF);
	fclose(f);
	char buf[0x100];
	sprintf(buf, "cpu %s &", fname);
	return system(buf);
}
