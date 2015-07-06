#include <stdio.h>
//#include <stdlib.h>

#include "well1024.c"

int main()
{
	well1024_seed(0);
	for (int i = 0; i < 1000; i++)
	{
		double r = well1024();
		printf("%g\n", r);
	}
	return 0;
}
