#ifndef XFONTS_ALL_H
#define XFONTS_ALL_H
#include "xfont_4x6.c"
#include "xfont_5x7.c"
#include "xfont_5x8.c"
#include "xfont_6x10.c"
#include "xfont_6x12.c"
#include "xfont_6x13.c"
#include "xfont_6x13B.c"
#include "xfont_6x13O.c"
#include "xfont_6x9.c"
#include "xfont_7x13.c"
#include "xfont_7x13B.c"
#include "xfont_7x13O.c"
//#include "xfont_7x14.c"
#include "xfont_7x14B.c"
#include "xfont_8x13.c"
#include "xfont_8x13B.c"
#include "xfont_8x13O.c"
#include "xfont_9x15.c"
#include "xfont_9x15B.c"
#include "xfont_9x18.c"
#include "xfont_9x18B.c"
#include "xfont_10x20.c"
//#include "xfont_12x13ja.c"
//#include "xfont_18x18ja.c"
//#include "xfont_18x18ko.c"
#include "xfont_canny.c"
//#include "xfont_cannyM.c"
#include "xfont_clR6x12.c"
#include "xfont_helvR12.c"

struct bitmap_font *xfonts_all[] = {
	xfont_4x6, xfont_5x7, xfont_5x8, xfont_6x10, xfont_6x12, xfont_6x13,
	xfont_6x13B, xfont_6x13O, xfont_6x9, xfont_7x13, xfont_7x13B,
	xfont_7x13O, xfont_7x14B, xfont_8x13, xfont_8x13B, xfont_8x13O,
	xfont_9x15, xfont_9x15B, xfont_9x18, xfont_9x18B, xfont_10x20,
	xfont_canny, xfont_clR6x12, xfont_helvR12, 0
};
#endif//XFONTS_ALL_H
