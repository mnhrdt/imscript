//Copyright (c) 2011 ashelly.myopenid.com under <http://www.opensource.org/licenses/mit-license>
//
// Image filtering adaptation by Gabriele Facciolo, 2017
//

#include <stdlib.h>
#include <math.h>

//Customize for your data Item type
typedef float Item;
#define ItemLess(a,b)  ((a)<(b))
#define ItemMean(a,b)  (((a)+(b))/2)

struct Mediator
{
   Item* data;  //circular queue of values
   int*  pos;   //index into `heap` for each value
   int*  heap;  //max/median/min heap holding indexes into `data`.
   int   N;     //allocated size.
   int   idx;   //position in circular queue
   int   ct;    //count of items in queue
};

/*--- Helper Functions ---*/

#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2)   //count of items in maxheap

//returns 1 if heap[i] < heap[j]
int mmless(struct Mediator* m, int i, int j)
{
   return ItemLess(m->data[m->heap[i]],m->data[m->heap[j]]);
}

//swaps items i&j in heap, maintains indexes
int mmexchange(struct Mediator* m, int i, int j)
{
   int t = m->heap[i];
   m->heap[i]=m->heap[j];
   m->heap[j]=t;
   m->pos[m->heap[i]]=i;
   m->pos[m->heap[j]]=j;
   return 1;
}

//swaps items i&j if i<j;  returns true if swapped
int mmCmpExch(struct Mediator* m, int i, int j)
{
   return (mmless(m,i,j) && mmexchange(m,i,j));
}

//maintains minheap property for all items below i/2.
void minSortDown(struct Mediator* m, int i)
{
   for (; i <= minCt(m); i*=2)
   {  if (i>1 && i < minCt(m) && mmless(m, i+1, i)) { ++i; }
      if (!mmCmpExch(m,i,i/2)) { break; }
   }
}

//maintains maxheap property for all items below i/2. (negative indexes)
void maxSortDown(struct Mediator* m, int i)
{
   for (; i >= -maxCt(m); i*=2)
   {  if (i<-1 && i > -maxCt(m) && mmless(m, i, i-1)) { --i; }
      if (!mmCmpExch(m,i/2,i)) { break; }
   }
}

//maintains minheap property for all items above i, including median
//returns true if median changed
int minSortUp(struct Mediator* m, int i)
{
   while (i>0 && mmCmpExch(m,i,i/2)) i/=2;
   return (i==0);
}

//maintains maxheap property for all items above i, including median
//returns true if median changed
int maxSortUp(struct Mediator* m, int i)
{
   while (i<0 && mmCmpExch(m,i/2,i))  i/=2;
   return (i==0);
}

/*--- Public Interface ---*/


//creates new Mediator: to calculate `nItems` running median.
//mallocs single block of memory, caller must free.
struct Mediator* MediatorNew(int nItems)
{
   int size = sizeof(struct Mediator)+nItems*(sizeof(Item)+sizeof(int)*2);
   struct Mediator* m=  malloc(size);
   m->data= (Item*)(m+1);
   m->pos = (int*) (m->data+nItems);
   m->heap = m->pos+nItems + (nItems/2); //points to middle of storage.
   m->N=nItems;
   m->ct = m->idx = 0;
   while (nItems--)  //set up initial heap fill pattern: median,max,min,max,...
   {  m->pos[nItems]= ((nItems+1)/2) * ((nItems&1)?-1:1);
      m->heap[m->pos[nItems]]=nItems;
   }
   return m;
}


//Inserts item, maintains median in O(lg nItems)
void MediatorInsert(struct Mediator* m, Item v)
{
   int isNew=(m->ct<m->N);
   int p = m->pos[m->idx];
   Item old = m->data[m->idx];
   m->data[m->idx]=v;
   m->idx = (m->idx+1) % m->N;
   m->ct+=isNew;
   if (p>0)         //new item is in minHeap
   {  if (!isNew && ItemLess(old,v)) { minSortDown(m,p*2);  }
      else if (minSortUp(m,p)) { maxSortDown(m,-1); }
   }
   else if (p<0)   //new item is in maxheap
   {  if (!isNew && ItemLess(v,old)) { maxSortDown(m,p*2); }
      else if (maxSortUp(m,p)) { minSortDown(m, 1); }
   }
   else            //new item is at median
   {  if (maxCt(m)) { maxSortDown(m,-1); }
      if (minCt(m)) { minSortDown(m, 1); }
   }
}

//returns median item (or average of 2 when item count is even)
Item MediatorMedian(struct Mediator* m)
{
   Item v= m->data[m->heap[0]];
   if ((m->ct&1)==0) { v= ItemMean(v,m->data[m->heap[-1]]); }
   return v;
}


/*--- Test Code ---*/
#include <stdio.h>
#include "iio.h"
//void PrintMaxHeap(struct Mediator* m)
//{
//   int i;
//   if(maxCt(m))
//      printf("Max: %3d",m->data[m->heap[-1]]);
//   for (i=2;i<=maxCt(m);++i)
//   {
//      printf("|%3d ",m->data[m->heap[-i]]);
//      if(++i<=maxCt(m)) printf("%3d",m->data[m->heap[-i]]);
//   }
//   printf("\n");
//}
//void PrintMinHeap(struct Mediator* m)
//{
//   int i;
//   if(minCt(m))
//      printf("Min: %3d",m->data[m->heap[1]]);
//   for (i=2;i<=minCt(m);++i)
//   {
//      printf("|%3d ",m->data[m->heap[i]]);
//      if(++i<=minCt(m)) printf("%3d",m->data[m->heap[i]]);
//   }
//   printf("\n");
//}
//
//void ShowTree(struct Mediator* m)
//{
//   PrintMaxHeap(m);
//   printf("Mid: %3d\n",m->data[m->heap[0]]);
//   PrintMinHeap(m);
//   printf("\n");
//}




// the type of a "getpixel" function
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel_1(float *I, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return I[i+j*w];
}

int main(int argc, char* argv[])
{
   if (argc<2) {
      fprintf(stderr, "fast median filtering using dxd windows:\n");
      fprintf(stderr, "usage: %s d [in [out]]\n", argv[0]);
      //                      0  1  2   3
      return 0;
   }
   char *filename_in  = argc > 2 ? argv[2] : "-";
   char *filename_out = argc > 3 ? argv[3] : "-";
   int size = atoi(argv[1]);

   int odd_size = (size/2)*2 + 1;
   int hsize = (odd_size-1)/2;
   if (hsize<1) {
      printf("filter %d is too small\n", hsize);
      return 0;
   }

   int w, h, pd;
   float *in = iio_read_image_float_split(filename_in, &w, &h, &pd);
   float *out = malloc(sizeof*out*w*h*pd);

   struct Mediator* m = MediatorNew(size*size);

   getpixel_operator p = getpixel_1;
   for (int c=0;c<pd;c++) {
      for(int j=0;j<h;j++) {
         for(int i=-hsize;i<w+hsize;i++) {
            for(int k=-hsize;k<=hsize;k++) {
               MediatorInsert(m, p(in + c*w*h, w, h, i, (j+k)) );
            }
            if((i-hsize)>=0 && (i-hsize)<w) {
               out[c*w*h +(j)*w +(i-hsize)]=MediatorMedian(m);
            }
         }
      }
   }

   iio_write_image_float_split(filename_out, out, w, h, pd);
   free(m);
   free(out);
   free(in);
}
