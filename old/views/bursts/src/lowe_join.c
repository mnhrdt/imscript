#include <stdio.h>
#include <stdlib.h>

int main()
{
	int next, curr = EOF;
	while ((next = getchar()) != EOF && next != '\n')
		;
	while ((next = getchar()) != EOF) {
		if (curr != EOF && (curr != '\n' || next != ' '))
			putchar(curr);
		curr = next;
	}
	putchar('\n');
	return EXIT_SUCCESS;
}
