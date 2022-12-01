#include <ctype.h>    // isprint
#include <stdio.h>    // getchar, printf
#include <termios.h>  // tcgetattr, tcsetattr

int main()
{
	// canonicalize terminal
	struct termios m, o;
	tcgetattr(0, &m);
	o = m;
	m.c_lflag &= !(ECHO | ICANON);
	tcsetattr(0, TCSANOW, &m);

	while (1)
	{
		fprintf(stderr, " (going to wait for a char)\n");
		int c = getchar();
		fprintf(stderr, " (got c=%d)\n", c);
		printf("got c = %d", c);
		if (isprint(c))
			printf(" '%c'", c);
		printf("\n");
		if (c == 27)
			break;
	}

	// un-canonicalize back
	tcsetattr(0, TCSANOW, &o);
	return 42;
}
