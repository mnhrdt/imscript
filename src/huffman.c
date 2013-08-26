#include <stdint.h>

// 1. API
//int huffman_encode(uint8_t *out, uint8_t *in, int n);
//int huffman_decode(uint8_t *out, uint8_t *in, int n);

// precondition "out" must contain enough pre-allocated space
// at worse, 260+n


// 2. IMPLEMENTATION


// a node in the huffman tree
struct hnode {
	int word;
	int freq;
	int mother, left, right;
};

#define HEAP_ENERGY(h,i) h[i].freq
#define HEAP_SWAP(h,i,j) do{struct hnode tmp_=h[i];h[i]=h[j];h[j]=tmp_;}while(0)
#include "abstract_heap.h"

int huffman_encode(uint8_t *out, uint8_t *in, int n, int nw)
{
#define NW nw
#define NT (NW*2-1)
	struct hnode queue[NW]; // the priority queue (starts with the leafs)
	struct hnode htree[NT]; // the Huffman tree (starts empty)

	// 0. initialize tree with empty data
	//for (int i = 0; i < NW; i++) {
	//	struct hnode *x = queue + i;
	//	x->word = x->freq = x->mother = x->left = x->right = -1;
	//}
	for (int i = 0; i < NT; i++) {
		struct hnode *x = htree + i;
		x->word = -1;
		x->mother = x->left = x->right = -1;
	}

	// 1. compute histogram
	for (int i = 0; i < NW; i++) {
		queue[i].word = i;
		queue[i].freq = 0;
		//queue[i].mother = -1;
	}
	for (int i = 0; i < n; i++) {
		assert(in[i] < NW);
		queue[in[i]].freq += 1;
	}

	HEAP_BUILD(queue, NW);

	for (int i = 0; i < NW - 1; i++)
	{
		struct hnode *qa = queue + NW - i - 1;
		HEAP_REMOVE_TOP(queue, 1 + qa - queue);
		//fprintf(stderr,"got top w=%d, freq=%d\n", qa->word, qa->freq);
		struct hnode *ta = htree + qa->word;

		struct hnode *qb = queue + NW - i - 2;
		HEAP_REMOVE_TOP(queue, 1 + qb - queue);
		//fprintf(stderr,"got top w=%d, freq=%d\n", qb->word, qb->freq);
		struct hnode *tb = htree + qb->word;

		struct hnode *tc = htree + NW + i;
		tc->left = qa->word;
		tc->right = qb->word;
		ta->mother = tb->mother = NW + i;

		qb->word = NW + i;
		qb->freq += qa->freq;
		//fprintf(stderr,"created w=%d, freq=%d\n", qb->word, qb->freq);
		HEAP_ADD(queue, qb - queue);

		//fprintf(stderr, "\n");
	}

	htree[NT-1].mother = NT-1;

	for (int i = 0; i < NT; i++)
	{
		struct hnode *x = htree + i;
		fprintf(stderr, "htree[%d] = {m=%d l=%d r=%d}\n",
			i, x->mother, x->left, x->right);
	}
	fprintf(stderr, "\n");


	return 0;
}

int main(void)
{
	uint8_t t[40] = {
		0, 1, 3,  3, 3, 3,  3, 3, 3,  0,
		2, 2, 2,  2, 2, 3,  3, 3, 3,  0,
		4, 4, 5,  5, 5, 6,  7, 7, 7,  7,
		7, 7, 7,  7, 7, 7,  6, 6, 6,  6,
	};
	uint8_t tt[100];
	return huffman_encode(tt, t, 40, 8);
}
