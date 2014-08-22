#include <stdio.h>

static void swap(char *t[], int i, int j)
{
	if (i != j)
	{
		char *tmp = t[i];
		t[i] = t[j];
		t[j] = tmp;
	}
}

static void emit(char *t[], int n)
{
	for (int i = 0; i < n; i++)
		printf("%s%c", t[i], i == n-1 ? '\n' : ' ');
}

static void emit_perms(char *t[], int s, int n)
{
	if (n > 0) {
		for (int i = 0; i < n; i++)
		{
			swap(t, s, s + i);
			emit_perms(t, s + 1, n - 1);
			swap(t, s, s + i);
		}
	} else
		emit(t, s);
}

int main(int c, char *v[])
{
	int n = c - 1;
	emit_perms(v + 1, 0, n);
	return 0;
}
