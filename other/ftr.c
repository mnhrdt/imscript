#define FTR_BACKEND_X11      1
#define FTR_BACKEND_WINDOWS  2
#define FTR_BACKEND_COCOA    3
#define FTR_BACKEND_GLUT     4
#define FTR_BACKEND_FREEGLUT 5
#define FTR_BACKEND_TEXT     6

#ifndef FTR_BACKEND
#  define FTR_BACKEND FTR_BACKEND_X11
#endif

#ifdef FTR_BACKEND

#  if FTR_BACKEND == FTR_BACKEND_X11
#    include "ftr_x11.c"
#  endif

#   if FTR_BACKEND == FTR_BACKEND_FREEGLUT
#    include "ftr_freeglut.c"
#   endif

#else

#   error "please define FTR BACKEND"

#endif
