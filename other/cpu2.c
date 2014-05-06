// gcc -std=c99 -O3 cpu2.c iio.o -lGL -lglut -lpng -ljpeg -ltiff

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <GL/glut.h>


//
// global state
//
//float *image;
//int image_w, image_h, image_pd;

struct state {
	// whole image
	int w, h, pd;
	float *image;

	// viewport
	int px, py, pw, ph, max_pw, max_ph;
	float *view;
//	unsigned char *view;
	float qmin, qmax;
	//float a, b;
};

struct state global_state[1];


// set-up functions
void set_up_global_state(void);
void set_up_opengl_environment(void);
void set_up_glut_environment(void);

// handlers
void handle_display(void);
void handle_reshape(int,int);
void handle_keyboard(unsigned char,int,int);
void handle_keyboard_special(int,int,int);
void handle_mouse(int,int,int,int);
void handle_motion_active(int,int);
void handle_motion_passive(int,int);

// actions
void action_color_lighten(void);
void action_color_darken(void);
void action_color_lighten_at(int, int);
void action_color_darken_at(int, int);
void action_color_center_at(int, int);
void action_viewport_translation(int, int);
void dump_view(void);


int main(int c, char **v)
{
	if (c < 2)
		return fprintf(stderr, "usage:\n\t%s image\n", *v);

	// read image into global variable
	float *iio_read_image_float_vec(char *, int*, int*, int*);
	int w, h, pd;
	float *x = iio_read_image_float_vec(v[1], &w, &h, &pd);
	if (pd != 3 && pd != 4)
		return fprintf(stderr, "only RGB or RGBA images!\n");
	global_state->w = w;
	global_state->h = h;
	global_state->pd = pd;
	global_state->image = x;


	// set up environment
	glutInit(&c, v);
	set_up_global_state();
	set_up_glut_environment();
	set_up_opengl_environment();

	// run
	glutMainLoop();

	// not reached
	return 42;
}


void set_up_global_state(void)
{
	struct state *e = global_state;

	e->max_pw = 1024;
	e->max_ph = 1024;

	e->pw = e->w < e->max_pw ? e->w : e->max_pw;
	e->ph = e->h < e->max_ph ? e->h : e->max_ph;

	e->px = 0;
	e->py = 0;
	//e->a = 1;
	//e->b = 0;
	e->qmin = 0;
	e->qmax = 255;
	e->view = malloc(e->max_pw * e->max_ph * 3 * sizeof*e->view);

	dump_view();
}

void set_up_opengl_environment(void)
{
	//glPixelTransferf(GL_RED_SCALE, 1.0/255);
	//glPixelTransferf(GL_GREEN_SCALE, 1.0/255);
	//glPixelTransferf(GL_BLUE_SCALE, 1.0/255);
}

void set_up_glut_environment(void)
{
	// set display mode
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	// set window
	glutCreateWindow("Whatever");
	glutReshapeWindow(global_state->pw, global_state->ph);

	// set window handlers
	glutDisplayFunc(handle_display);
	glutReshapeFunc(handle_reshape);
	glutIdleFunc(NULL);

	// set mouse handlers
	glutMouseFunc(handle_mouse);
	glutMotionFunc(handle_motion_active);
	glutPassiveMotionFunc(handle_motion_passive);

	// set keyboard handlers
	glutKeyboardFunc(handle_keyboard);
	glutSpecialFunc(handle_keyboard_special);
}


// event: resize window
void handle_reshape(int pw, int ph)
{
	struct state *e = global_state;

	if (pw > e->max_pw) pw = e->max_pw;
	if (ph > e->max_ph) ph = e->max_ph;
	e->pw = pw;
	e->ph = ph;
	glViewport(0, 0, pw, ph);

	dump_view();
	glutPostRedisplay();
}

// event: redraw window
void handle_display(void)
{
	glRasterPos2i(-1,1);
	glPixelZoom(1,-1);

	int w = global_state->pw;
	int h = global_state->ph;
	//unsigned char *x = global_state->view;
	//glDrawPixels(w, h, GL_RGB, GL_UNSIGNED_BYTE, x);
	glDrawPixels(w, h, GL_RGB, GL_FLOAT, global_state->view);

	glutSwapBuffers();
}

// event: regular keypress
void handle_keyboard(unsigned char key,int x, int y)
{
	fprintf(stderr, "key %d '%c' (%d %d)\n", key, key, x, y);
	switch (key) {
	case 27: // ESC
	case 'Q':
	case 'q': 
		exit(0); 
		break;
	}
}

// event: special key (arrow, function, etc)
void handle_keyboard_special(int key, int x, int y)
{
	//fprintf(stderr, "special key %d (%d %d)\n", key, x, y);

	int d = 1 * (1 << (glutGetModifiers()*2));

	switch (key) {
	case GLUT_KEY_UP: action_viewport_translation(0, d); break;
	case GLUT_KEY_DOWN: action_viewport_translation(0, -d); break;
	case GLUT_KEY_LEFT: action_viewport_translation(-d, 0); break;
	case GLUT_KEY_RIGHT: action_viewport_translation(d, 0); break;
	default: return;
	}

	dump_view();
	glutPostRedisplay();
}

// event: mouse CLICK
void handle_mouse(int button, int state,int x,int y)
{
	// button=0 : left
	// button=1 : middle
	// button=2 : right (not seen becauset it is captured by menu)
	// button=3 : wheel up
	// button=4 : wheel down

	if (state == GLUT_DOWN) {
		int m = glutGetModifiers(); // 1=shift, 2=ctrl, 4=alt
		fprintf(stderr, "mouse but=%d, state=%d, pos=%d %d mod=%d\n",
				button, state, x, y, m);

		if (m == 0 && button == 3) action_color_lighten();
		if (m == 0 && button == 4) action_color_darken();
		if (m == 1 && button == 3) action_color_lighten_at(x, y);
		if (m == 1 && button == 4) action_color_darken_at(x, y);

		dump_view();
		glutPostRedisplay();
	}
}

// event: mouse DRAG
void handle_motion_active(int x,int y)
{
	int m = glutGetModifiers(); // 1=shift, 2=ctrl, 4=alt
	fprintf(stderr, "motion active %d,%d (m=%d)\n", x, y, m);
}

// event: mouse MOTION
void handle_motion_passive(int x, int y)
{
	int m = glutGetModifiers(); // 1=shift, 2=ctrl, 4=alt

	if (m == 1) {
		action_color_center_at(x, y);

		dump_view();
		glutPostRedisplay();
	}
}



void action_color_lighten(void)
{
	global_state->qmin /= 1.1;
	global_state->qmax /= 1.1;
}

void action_color_darken(void)
{
	//global_state->a /= 1.1;
	global_state->qmin *= 1.1;
	global_state->qmax *= 1.1;
}

float get_avg_below(int x, int y)
{
	struct state *e = global_state;
	x += e->px;
	y += e->py;
	if (x < 0 || y < 0 || x >= e->w || y >= e->h)
		return 0;
	float avg = 0;
	for (int l = 0; l < 3; l++)
		avg += e->image[(e->w*y+x)*e->pd+l];
	return avg/3;
}

void action_color_lighten_at(int x, int y)
{
	struct state *e = global_state;
	float avg = get_avg_below(x, y);
	fprintf(stderr, "stderr action color LIGHTEN at %d %d (%g)\n",x,y,avg);

	float hrange = (e->qmax - e->qmin)/2;
	float center = (e->qmax + e->qmin)/2;
	hrange *= 1.1;
	e->qmin = center - hrange;
	e->qmax = center + hrange;
}

void action_color_darken_at(int x, int y)
{
	struct state *e = global_state;
	float avg = get_avg_below(x, y);
	fprintf(stderr, "stderr action color DARKEN at %d %d (%g)\n",x,y,avg);

	float hrange = (e->qmax - e->qmin)/2;
	float center = (e->qmax + e->qmin)/2;
	hrange /= 1.1;
	e->qmin = center - hrange;
	e->qmax = center + hrange;
}

void action_color_center_at(int x, int y)
{
	struct state *e = global_state;
	float avg = get_avg_below(x, y);
	//fprintf(stderr, "stderr action color CENTER at %d %d (%g)\n",x,y,avg);

	float hrange = (e->qmax - e->qmin)/2;
	e->qmin = avg - hrange;
	e->qmax = avg + hrange;
}

void action_viewport_translation(int dx, int dy)
{
	global_state->px += dx;
	global_state->py -= dy;
}

void dump_view(void)
{
	struct state *e = global_state;
	memset(e->view, 0, 3 * e->pw * e->ph * sizeof*e->view);
	for (int j = 0; j < e->ph; j++)
	for (int i = 0; i < e->pw; i++)
	for (int l = 0; l < 3; l++)
	{
		int ii = i + e->px;
		int jj = j + e->py;
		if (ii < 0 || ii >= e->w) continue;
		if (jj < 0 || jj >= e->h) continue;
		int idx = (j*e->pw + i)*3 + l;
		int iidx = (jj*e->w + ii)*e->pd + l;
		float v = e->image[iidx];
		e->view[idx] = (v - e->qmin)/(e->qmax - e->qmin);
		//unsigned char r;
		//if (v < e->qmin) r = 0;
		//else if (v > e->qmax) r = 255;
		//else r = 255*(v - e->qmin)/(e->qmax - e->qmin);
		//e->view[idx] = r;
		//
		//float v = e->a * e->image[iidx] + e->b;
		//e->view[idx] = fmax(0,fmin(255,v));
	}
	//fprintf(stderr, "DUMP (%d %d %d %d) qmin qmax = %g %g\n",
	//		e->px, e->py, e->pw, e->ph,
	//		e->qmin, e->qmax);
}
