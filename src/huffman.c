#include <stdint.h>

// 1. API
//int huffman_encode(uint8_t *out, uint8_t *in, int n);
//int huffman_decode(uint8_t *out, uint8_t *in, int n);

// precondition "out" must contain enough pre-allocated space
// at worse, 260+n


// 2. IMPLEMENTATION
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


struct tree_node { int mother, left, right; };
struct heap_node { int word, freq; };

struct code_word { int word, length, *bits; };

static void fill_code_word(struct code_word *w, struct tree_node *t, int l)
{
	w->word = l;
	int len = 0, i = l, cx = 1;
	while (t[i].mother != i) {
		len += 1;
		i = t[i].mother;
	}
	w->length = len;
	w->bits = malloc((w->length+1) * sizeof*w->bits);
	i = l;
	while (t[i].mother != i) {
		int m = t[i].mother;
		//int bit = (i == t[m].left) ? 1 : 0;
		//int bitpos = w->length - cx++;
		//fprintf(stderr, "word %d, bit %d = %d\n", l, bitpos, bit);
		//w->bits[bitpos] = bit;
		w->bits[w->length-cx++] = (i == t[m].left) ? 1 : 0;
		i = m;
	}
}

//static void print_code_of_this_leaf(struct tree_node *t, int n, int l)
//{
//	assert(l >= 0);
//	assert(l < n);
//	int len = 0, i = l;
//	while (t[i].mother != i)
//	{
//		len += 1;
//		i = t[i].mother;
//	}
//	//fprintf(stderr, "len(%d) = %d\n", l, len+1);
//	int code[len], cx = 0;
//	i = l;
//	while (t[i].mother != i)
//	{
//		int m = t[i].mother;
//		code[len-++cx] = (i == t[m].left) ? 1 : 0;
//		i = m;
//	}
//	//fprintf(stderr, "cx = %d\n", cx);
//	fprintf(stderr, "code(%d) = 1", l);
//	for (int j = 0; j < len; j++)
//		fprintf(stderr, " %d", code[j]);
//	fprintf(stderr, "\n");
//
//}

static void find_code_from_tree(struct code_word *c, struct tree_node *t, int n)
{
	assert(1 == n%2);
	int nw = (n + 1)/2;
	for (int i = 0; i < nw; i++)
	{
		assert(t[i].left < 0 || t[i].left == i);
		fill_code_word(c + i, t, i);
	}
	for (int i = 0; i < nw; i++)
	{
		struct code_word *w = c + i;
		fprintf(stderr, "codeword %d:", w->word);
		for (int j = 0; j < w->length; j++)
			fprintf(stderr, " %d", w->bits[j]);
		fprintf(stderr, "\n");
	}
}

#define HEAP_ENERGY(h,i) h[i].freq
#define HEAP_SWAP(h,i,j) do{struct heap_node t_=h[i];h[i]=h[j];h[j]=t_;}while(0)
#include "abstract_heap.h"

static
void compute_huffman_tree(struct tree_node *htree, uint8_t *in, int n, int nw)
{
#define NW nw
#define NT (NW*2-1)
	struct heap_node queue[NW]; // the priority queue

	// 0. initialize the tree
	for (int i = 0; i < NT; i++)
		htree[i].mother = htree[i].left = htree[i].right = i;

	// 1. compute histogram
	for (int i = 0; i < NW; i++) {
		queue[i].word = i;
		queue[i].freq = 0;
	}
	for (int i = 0; i < n; i++) {
		assert(in[i] < NW);
		queue[in[i]].freq += 1;
	}

	HEAP_BUILD(queue, NW);

	for (int i = 0; i < NW - 1; i++)
	{
		struct heap_node *qa = queue + NW - i - 1;
		HEAP_REMOVE_TOP(queue, 1 + qa - queue);
		//fprintf(stderr,"got top w=%d, freq=%d\n", qa->word, qa->freq);
		struct tree_node *ta = htree + qa->word;

		struct heap_node *qb = queue + NW - i - 2;
		HEAP_REMOVE_TOP(queue, 1 + qb - queue);
		//fprintf(stderr,"got top w=%d, freq=%d\n", qb->word, qb->freq);
		struct tree_node *tb = htree + qb->word;

		struct tree_node *tc = htree + NW + i;
		tc->left = qa->word;
		tc->right = qb->word;
		ta->mother = tb->mother = NW + i;

		qb->word = NW + i;
		qb->freq += qa->freq;
		//fprintf(stderr,"created w=%d, freq=%d\n", qb->word, qb->freq);
		HEAP_ADD(queue, qb - queue);

		//fprintf(stderr, "\n");
	}

	//htree[NT-1].mother = NT-1;

	//for (int i = 0; i < NT; i++)
	//{
	//	struct tree_node *x = htree + i;
	//	fprintf(stderr, "htree[%d] = {m=%d l=%d r=%d}\n",
	//		i, x->mother, x->left, x->right);
	//}
	//fprintf(stderr, "\n");
}

#define SETBIT(x,i) ((x)|=(1<<(i)))
#define GETBIT(x,i) (bool)((x)&(1<<(i)))

static void add_bit_to_stream(uint8_t *x, int n, bool value)
{
	assert(value == 0 || value == 1);
	int pos = n/8;
	int res = n%8;
	if (res == 1)
		x[pos] = 0;
	if (value)
		SETBIT(x[pos],res);
}



static bool get_bit_from_stream(uint8_t *x, int n)
{
	int pos = n/8;
	int res = n%8;
	return GETBIT(x[pos], res);
}

static
int encode_message(uint8_t *y, struct code_word *c, int nw, uint8_t *x, int n)
{
	assert(nw < 256);
	fprintf(stderr, "original message of %d bytes\n", n);

	// 1. encode the codebook (each int goes to 4 bytes)
	uint32_t *yy = (void*)y;
	// 1.1. number of words
	yy[0] = nw;
	// 1.2. the lengths of each word
	for (int i = 0; i < nw; i++)
		yy[i] = c[i].length;
	int r = 4*nw + 4, nbits = 0;
	y += 4*(nw + 1);
	// 1.3. the words themselves
	for (int i = 0; i < nw; i++)
		for (int j = 0; j < c[i].length; j++)
			add_bit_to_stream(y, nbits++, c[i].bits[j]);
	fprintf(stderr, "header of %d bytes\n", r);

	// 2. encode the message
	for (int i = 0; i < n; i++)
		for (int j = 0; j < c[x[i]].length; j++)
			add_bit_to_stream(y, nbits++, c[x[i]].bits[j]);

	int nbytes = nbits/8 + (bool)(nbits%8);
	r += nbytes;
	fprintf(stderr, "message of %d bytes\n", nbytes);
	fprintf(stderr, "for a total of %d bytes (%g%%)\n", r, 100.0*r/n);
	return r;
}

int huffman_encode(uint8_t *out, uint8_t *in, int n, int nw)
{
	struct tree_node htree[NT];
	compute_huffman_tree(htree, in, n, nw);

	struct code_word code[NW];
	find_code_from_tree(code, htree, NT);

	return encode_message(out, code, NW, in, n);
}

int main(int c, char **v)
{
	if (c != 3) return -1;
	int n = atoi(v[1]);
	int m = atoi(v[2]);
	uint8_t t[n], tt[2*n+1000];
	for (int i = 0; i < n; i++)
		t[i] = fmin(rand()%m,rand()%m);//sqrt(1.0*(rand()%m)/m)*m);
	return huffman_encode(tt, t, n, m);
	//uint8_t t[40] = {
	//	0, 1, 3,  3, 3, 3,  3, 3, 3,  0,
	//	2, 2, 2,  2, 2, 3,  3, 3, 3,  0,
	//	4, 4, 5,  5, 5, 6,  7, 7, 7,  7,
	//	7, 7, 7,  7, 7, 7,  6, 6, 6,  6,
	//};
	//uint8_t tt[100];
	//return huffman_encode(tt, t, 40, 8);
}
