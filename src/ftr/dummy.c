#include <stdio.h>
#include <termios.h>

int main()
{
	struct termios m;
	tcgetattr(0, &m);
	m.c_lflag &= !(ECHO | ICANON);
	tcsetattr(0, TCSANOW, &m);

	while (1)
	{
		int c = getchar();
		printf("got c = %d\n", c);
		if (c == 27)
			break;
	}

	return 0;
}
