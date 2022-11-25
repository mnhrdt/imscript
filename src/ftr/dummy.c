#include <stdio.h>
#include <termios.h>

int main()
{
	struct termios m, o;
	tcgetattr(0, &m);
	o = m;
	m.c_lflag &= !(ECHO | ICANON);
	tcsetattr(0, TCSANOW, &m);

	while (1)
	{
		int c = getchar();
		printf("got c = %d\n", c);
		if (c == 27)
			break;
	}

	tcsetattr(0, TCSANOW, &o);
	return 42;
}
