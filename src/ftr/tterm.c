// some trials for tty-based event loop

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <termios.h>

static int gettsize(int *w, int *h)
{
	char *b;
	b = getenv("COLUMNS"); if (!b) return 0; *w = atoi(b);
	b = getenv("LINES");   if (!b) return 0; *h = atoi(b);
	return 2;
}

int main(void)
{
	// disable canonical and echo modes
	struct termios t[1];
	tcgetattr(0, t);
	t->c_lflag &= ~ICANON;
	t->c_lflag &= ~ECHO;
	tcsetattr(0, TCSANOW, t);

	// event loop
	int buf_n = 0x100, buf[buf_n], buf_i = 0;
	while (1)
	{
		int c = getchar();
		if (c == EOF) break;

		int w, h;
		int r = gettsize(&w, &h);
		if (!r) w = h = -1;
		printf("GOT{%dx%d} %d '%c'\n", w, h, c, isprint(c) ? c : ' ');
		if (c == 'q')
			break;
	}

	// cleanup and exit
	tcgetattr(0, t);
	t->c_lflag |= ICANON;
	t->c_lflag |= ECHO;
	tcsetattr(0, TCSANOW, t);
	return 0;
}
