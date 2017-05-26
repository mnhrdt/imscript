#ifndef _ABSTRACT_HEAP_H
#define _ABSTRACT_HEAP_H

// See a comment labelled "Documentation:" below



#ifndef assert
#include <assert.h>
#endif
#include <stdio.h>

// Note that HEAP_ENERGY can be a sequence of expressions, e.g., such as
// "assert(i>=0&&i<LENGTH(h)),h[i]".  It MUST NOT be enclosed on parentheses.
// (this is a weird hack on the precedence of the "," and "=" operators)
//
// 	fuck(h,i), h[i] = e;	// fucks, and then assigns "h[i] := e"
// 	e = (fuck(h,i), h[i]);	// fucks, and then assigns "e := h[i]"

#ifndef HEAP_ENERGY
//#define HEAP_ENERGY(h,i) h[i]
#error "you must define HEAP_ENERGY to use a heap!"
#endif

#ifndef HEAP_SWAP
#ifdef __GNUC__
#define HEAP_SWAP(h,i,j) do{__typeof__(*h) t_=h[i];h[i]=h[j];h[j]=t_;}while(0)
#else
#error "you must define HEAP_SWAP to use an abstract heap!"
#endif // __GNUC__
#endif // HEAP_SWAP

#ifndef HEAP_HIGH_ENERGY
#define CMPEQ <=
#define CMPIN <
#else
#define CMPEQ >=
#define CMPIN >
#endif // HEAP_HIGH_ENERGY


// Documentation:
//
// These macros act over an arbitrary given array.  Some macros are
// intended to be used internally only.  The "interface" looks like
// this (in the case that the array is a common array of floating-point
// values):
//
// #define HEAP_ENERGY(h,i) h[i]
// #define HEAP_SWAP(h,i,j) do{double t_t=h[i];h[i]=h[j];h[j]=t_t;}while(0)
// #include "abstract_heap.h"
//
// ...
// 	// build an array of random floats and sort it:
// 	int n = 100;
// 	float *t=malloc(n * sizeof(*t));
// 	for (int i = 0; i < n; i++)
// 		t[i] = rand();
// 	HEAP_SORT(t,n);
// 	// now the values of t are sorted (and t is not a heap anymore!)
//
//	// other available operations:
//	HEAP_BUILD(h,n);	// build a heap out of an arbitrary array
//	h[0];			// access the top of the heap
//	HEAP_REMOVE_TOP(h,n);	// put the top at the end and re-arrange
//	n -= 1;			// this is needed after removal of the top
//	HEAP_CHANGE_ENERGY(h,n,i,E);	// change the energy of ith element
//					// and re-arrange
//
//
// The main interest of this approach is that the elements of the array can be
// arbitrary data structures, and the user-defined function HEAP_SWAP can
// perform a lot of neat tricks, like keeping track of pointers that refer to
// the elements of the heap from outside, who get updated gracefully as the
// heap rearranges its elements.  Try that on C++! (no, overloading the copy
// operator to change its lvalue is NOT considered elegant).
//
//

// these three functions are only intended for internal use
#define HEAP_PARENT(x) (((x)-1)/2)
#define HEAP_LEFT(x) (2*(x)+1)
#define HEAP_RIGHT(x) (2*((x)+1))

// assert the heap condition
#define HEAP_ASSERT(h,n) do{\
	for(int _i_as = 1; _i_as < (n); _i_as++)\
	{\
		assert((HEAP_ENERGY(h,HEAP_PARENT(_i_as)))\
				CMPEQ (HEAP_ENERGY(h,_i_as)));\
	}\
}while(0)

// internal function
#define HEAP_FIXUP(h,n,i_fu) do{\
        int _iifu=i_fu;\
        while(_iifu && \
		(HEAP_ENERGY(h,_iifu)) CMPIN (HEAP_ENERGY(h,HEAP_PARENT(_iifu)))\
	     )\
        {\
                HEAP_SWAP(h,_iifu,HEAP_PARENT(_iifu));\
                _iifu=HEAP_PARENT(_iifu);\
        }\
}while(0)

// USER-VISIBLE FUNCTION: build a heap out of an arbitrary array
#define HEAP_BUILD(h,n) do{\
	int _i_bu = 0;\
	while(_i_bu < (n))\
	{\
		HEAP_FIXUP(h, (n), _i_bu);\
		_i_bu++;\
	}\
}while(0)

// internal function, (best-of-three)
#define HEAP_BOT(h,n,i_bo) (HEAP_LEFT(i_bo)>=(n))?(i_bo):(\
		(HEAP_RIGHT(i_bo)==(n))?((HEAP_ENERGY(h,i_bo))CMPIN (HEAP_ENERGY(h,HEAP_LEFT(i_bo)))?i_bo:HEAP_LEFT(i_bo)):(\
		((HEAP_ENERGY(h,i_bo))CMPEQ (HEAP_ENERGY(h,HEAP_LEFT(i_bo)))&&(HEAP_ENERGY(h,i_bo))CMPEQ (HEAP_ENERGY(h,HEAP_RIGHT(i_bo))))?(i_bo):(\
		((HEAP_ENERGY(h,HEAP_LEFT(i_bo)))CMPEQ (HEAP_ENERGY(h,i_bo))&&(HEAP_ENERGY(h,HEAP_LEFT(i_bo)))CMPEQ (HEAP_ENERGY(h,HEAP_RIGHT(i_bo))))?(HEAP_LEFT(i_bo)):(\
		((HEAP_ENERGY(h,HEAP_RIGHT(i_bo)))CMPEQ (HEAP_ENERGY(h,HEAP_LEFT(i_bo)))&&(HEAP_ENERGY(h,HEAP_RIGHT(i_bo)))CMPEQ (HEAP_ENERGY(h,i_bo)))?(HEAP_RIGHT(i_bo)):-1))))

// internal function
#define HEAP_FIXDOWN(h,n,i_fd) do{\
	int _j_fd;\
	int _i_fdl = i_fd;\
	while((_j_fd=HEAP_BOT(h,n,_i_fdl)) != _i_fdl)\
	{\
		HEAP_SWAP(h, _i_fdl, _j_fd);\
		_i_fdl = _j_fd;\
	}\
}while(0)

// USER-VISIBLE FUNCTION:
// (n has to be decreased manually by the user afterwards)
#define HEAP_REMOVE_TOP(h,n) do{\
	HEAP_SWAP(h, 0, (n) - 1);\
	HEAP_FIXDOWN(h, (n) - 1, 0);\
}while(0)

// USER-VISIBLE FUNCTION:
// the energy of the new element has to be set beforehand
// "k" is the location of the new element
// (so the heap gets to have "k+1" elements)
#define HEAP_ADD(h,k) do{\
	HEAP_FIXUP(h, k+1, k);\
}while(0)

// USER-VISIBLE FUNCTION: sorts an arbitrary array
#define HEAP_SORT(h,n) do{\
	int _n_so = n;\
	HEAP_BUILD(h, _n_so);\
	while(_n_so > 1)\
	{\
		HEAP_SWAP(h, 0, _n_so - 1);\
		HEAP_FIXDOWN(h, _n_so - 1, 0);\
		_n_so -= 1;\
	}\
}while(0)

// USER-VISIBLE FUNCTION: update the energy of an element
#define HEAP_CHANGE_ENERGY(h,n,i_ce,E) do{\
	double _oE_ce = (HEAP_ENERGY(h, i_ce));\
	HEAP_ENERGY(h, i_ce) = E;\
	if (_oE_ce CMPIN E)\
		HEAP_FIXDOWN(h, n, i_ce);\
	else\
		HEAP_FIXUP(h, n, i_ce);\
}while(0)

// display a small heap on stdout (assumes energies are doubles)
#define HEAP_DISPLAY(h,n) do{\
	printf("heap %d\n",(n));\
	for(int _i_di = 0; _i_di < (n); _i_di++)\
		printf("%d:\t%g\n", _i_di, (HEAP_ENERGY(h,_i_di)));\
}while(0)


#endif /* _ABSTRACT_HEAP_H */
