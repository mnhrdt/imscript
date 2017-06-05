// sort algorithms

typedef void (sort_function_ints)(int *, int);

#include <assert.h>


static void swap_places(int *x, int a, int b)
{
	if (a == b) return;
	//printf("SWAP %d %d\n", a, b);
	int t = x[a];
	x[a] = x[b];
	x[b] = t;
}

void quicksort_ints(int *t, int n)
{
	if (n < 2) return;
	int p = n/2;
	swap_places(t, 0, p);
	int top = 1;
	for (int i = 1; i < n; i++)
		if (t[i] < t[0])
			swap_places(t, i, top++);
	swap_places(t, 0, top-1);
	quicksort_ints(t, top-1);
	quicksort_ints(t+top, n-top);
}

void mergesort_ints(int *t, int n)
{
	if (n < 2) return;
	int *ta = t,     na = n/2,   ia = 0;
	int *tb = t+n/2, nb = n-n/2, ib = 0;
	mergesort_ints(ta, na);
	mergesort_ints(tb, nb);
	int tmp[n], cx = 0;
	while (ia < na && ib < nb)
		if (ta[ia] < tb[ib])
			tmp[cx++] = ta[ia++];
		else
			tmp[cx++] = tb[ib++];
	if (ia < na)
		while (ia < na)
			tmp[cx++] = ta[ia++];
	else
		while (ib < nb)
			tmp[cx++] = tb[ib++];
	for (int i = 0; i < n; i++)
		t[i] = tmp[i];
}



#define HEAP_PARENT(i) (((i)-1)/2)
#define HEAP_LEFT(i) (2*(i)+1)
#define HEAP_RIGHT(i) (2*(i)+2)

static void max_heap_fix_up(int *t, int p)
{
	while (p && t[p] > t[HEAP_PARENT(p)])
	{
		swap_places(t, p, HEAP_PARENT(p));
		p = HEAP_PARENT(p);
	}

}

static void max_heap_fix_down(int *t, int n, int p)
{
	while (1)
	{
		int i = p;
		int left = HEAP_LEFT(i);
		int right = HEAP_RIGHT(i);
		if (left < n && t[i] < t[left])
			i = left;
		if (right < n && t[i] < t[right])
			i = right;
		if (i == p)
			break;
		swap_places(t, p, i);
		p = i;
	}
}

void heapsort_ints(int *t, int n)
{
	// build max-heap
	for (int i = 1; i < n; i++)
		max_heap_fix_up(t, i);

	// move top elements to the end
	for (int i = 1; i < n; i++)
	{
		swap_places(t, 0, n-i);
		max_heap_fix_down(t, n-i, 0);
	}
}


//static void insert(int *t, int n, int x)
//{
//	// skip elements less than x
//	for (int i = 0; i
//}
//
//void insertsort_ints(int *t, int n)
//{
//	for (int i = 1; i < n; i++)
//		insert(t, i, t[i]);
//}

#include <stdio.h>
#include <stdlib.h>

int main(int c, char *v[])
{
	if (c != 3)
		return 1;
	int n = atoi(v[1]);
	int m = atoi(v[2]);
	int t[n];
	for (int i = 0; i < n; i++)
	{
		t[i] = rand()%m;
		printf("before[%d] = %d\n", i, t[i]);
	}
	//sort_function_ints *f = mergesort_ints;
	//sort_function_ints *f = quicksort_ints;
	sort_function_ints *f = heapsort_ints;
	f(t, n);
	for (int i = 0; i < n; i++)
		printf("after[%d] = %d\n", i, t[i]);
	return 0;
}
