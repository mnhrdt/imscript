// gcc -std=c99 -O3 cpu.c iio.o -lGL -lglut -lpng -ljpeg -ltiff

#include <stdio.h>
#include <GL/glut.h>


//
// global state
//
float *image;
int image_w, image_h, image_pd;


// set-up functions
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


int main(int c, char **v)
{
	if (c < 2)
		return fprintf(stderr, "usage:\n\t%s image\n", *v);

	// read image into global variable
	float *iio_read_image_float_vec(char *, int*, int*, int*);
	image = iio_read_image_float_vec(v[1], &image_w, &image_h, &image_pd);
	if (image_pd != 3)
		return fprintf(stderr, "only RGB images!\n");

	// set up environment
	glutInit(&c, v);
	set_up_glut_environment();
	set_up_opengl_environment();

	// run
	glutMainLoop();

	// not reached
	return 42;
}


void set_up_opengl_environment(void)
{
	glPixelTransferf(GL_RED_SCALE, 1.0/255);
	glPixelTransferf(GL_GREEN_SCALE, 1.0/255);
	glPixelTransferf(GL_BLUE_SCALE, 1.0/255);
}

void set_up_glut_environment(void)
{
	// set display mode
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	// set window
	glutCreateWindow("Whatever");
	glutReshapeWindow(image_w, image_h);

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
void handle_reshape(int w, int h)
{
	glViewport(0, 0, w, h);
}

// event: redraw window
void handle_display(void)
{
	glRasterPos2i(-1,1);
	glPixelZoom(1,-1);

	glDrawPixels(image_w, image_h, GL_RGB, GL_FLOAT, image);

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
void handle_keyboard_special(int key,int x, int y)
{
	fprintf(stderr, "special key %d (%d %d)\n", key, x, y);
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

		if (button == 4 && m==0) // make image lighter
			for (int j = 0; j < image_h; j++)
			for (int i = 0; i < image_w; i++)
			for (int l = 0; l < image_pd; l++)
			{
				int idx = (image_w*j + i)*image_pd + l;
				image[idx] *= 1.1;
			}
		if (button == 3 && m==0) // make image darker
			for (int j = 0; j < image_h; j++)
			for (int i = 0; i < image_w; i++)
			for (int l = 0; l < image_pd; l++)
			{
				int idx = (image_w*j + i)*image_pd + l;
				image[idx] /= 1.1;
			}
		if ((button==4 || button==3 ) && m==1) // change local contrast
		{
			float factor = button == 3 ? 1.0/1.1 : 1.1;
			if (x < 0 || y < 0 || x >= image_w || y >= image_h)
				return;
			float col[3];
			for (int l = 0; l < 3; l++)
				col[l] = image[(image_w*y+x)*image_pd+l];
			for (int j = 0; j < image_h; j++)
			for (int i = 0; i < image_w; i++)
			for (int l = 0; l < image_pd; l++)
			{
				int idx = (image_w*j + i)*image_pd + l;
				image[idx] = factor*(image[idx]-col[l])+col[l];
			}
		}

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
void handle_motion_passive(int x,int y)
{
	int m = glutGetModifiers(); // 1=shift, 2=ctrl, 4=alt
	if (m) {
		if (m == 1) {
			if (x < 0 || y < 0 || x >= image_w || y >= image_h)
				return;
			float col[3];
			for (int l = 0; l < 3; l++)
				col[l] = image[(image_w*y+x)*image_pd+l];
			for (int j = 0; j < image_h; j++)
			for (int i = 0; i < image_w; i++)
			for (int l = 0; l < image_pd; l++)
			{
				int idx = (image_w*j + i)*image_pd + l;
				image[idx] += 127 - col[l];
			}
			glutPostRedisplay();
		} 
	}
}
