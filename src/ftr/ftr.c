#define FTR_BACKEND_X11      'x'
#define FTR_BACKEND_WINDOWS  'w'
#define FTR_BACKEND_COCOA    'c'
#define FTR_BACKEND_GLUT     'g'
#define FTR_BACKEND_FREEGLUT 'f'
#define FTR_BACKEND_TEXT     't'


#ifndef FTR_BACKEND
#  define FTR_BACKEND FTR_BACKEND_X11
#endif

#if FTR_BACKEND == FTR_BACKEND_X11
#  include "ftr_x11.c"
#elif FTR_BACKEND == FTR_BACKEND_FREEGLUT
#  include "ftr_freeglut.c"
#else
#  error "please define FTR_BACKEND to a 'x' or 'f'"
#endif
