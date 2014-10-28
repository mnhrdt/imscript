#include <stdio.h>

static void swap(char *t[], int i, int j)
{
	char *tmp = t[i];
	t[i] = t[j];
	t[j] = tmp;
}

static void print_array(char *t[], int n)
{
	for (int i = 0; i < n; i++)
		printf("%s%c", t[i], i == n-1 ? '\n' : ' ');
}

static void perm(char *t[], int s, int n)
{
	for (int i = 0; i < n; i++)
	{
		swap(t, s    , s + i);
		perm(t, s + 1, n - 1);
		swap(t, s    , s + i);
	}
	if (n == 0)
		print_array(t, s);
}

int main(int c, char *v[])
{
	perm(v + 1, 0, c - 1);
	return 0;
}
