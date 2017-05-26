// This file uses Vim folds.  If you don't like folds, type "zR"

/* a vi-like program for 3D image visualization and editing {{{1
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *                     *************************
 *                     **  G L S H A P E S T  **
 *                     *************************
 *
 *               Tree of SHAPES visualizer using openGL
 *
 *        (a vi-like program for 3D image visualization and editing)
 *
 *
 *
 * usage: read a 3D tree of shapes from stdin (and optionally specify a
 * coloring image as an argument)
 * 	j k	child, parent shape
 * 	J K	bottom, top of branch (XXX: UNIMPLEMENTED)
 * 	u i	child, and parent (iterated 5 times)
 * 	h l	previous, next sibling (ordered by volume)
 * 	o	toggle surface orientation
 * 	s	toggle wireframe/solid
 * 	g G	increase/decrease pallette scaling
 * 	0 1 2	set marching cubes option (standard / random / smart)
 * 	m , .	Mumford-Shah colorization, increase, decrease lambda
 * 	q	quit program
 *
 *
 * see the "global_static_table_of_key_functions" for all the keybindings
 * see the "main()" function for the possible arguments and options
 *
 *
 *
 * wanted features:
 * 	print all the descriptors of the current shape
 * 	pick the maximal meaningful part of this shape (and the following ones)
 * 	DONE "select" a shape
 * 	DONE "select" a part of a shape
 * 	DONE save the current selection to a .coh file
 * 	UGLY show slices (and typical projections) of data
 * 	choose a nice and meaningful pallette (and legend!)
 * 	indicate current mode (for quiche-eaters)
 * 	...
 *
 * TODO:
 * 	document the modes
 * 	unspaghettify the code
 * 		DONE do not scale the surface, use a transformation matrix
 * 		PARTIAL separate colorization from interface
 * 	design sane interface
 * 	rationalize the usage matrix transforms (now it's ad-hoc an inflexible)
 * 	cache the call lists themselves, not the raw triangulations (maybe?)
 * 		(no!: they may need to be mumford-shahized afterwards)
 * 		((answer: cache both))
 * 	general plane-drawing code (per fer destrals!!!)
 *
 * NOTE: the colorization image given as input must have the rank [0:255]
 *
 *
 */

/* #includes {{{1 */
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <locale.h>
#include <getopt.h>


#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#include <gtk/gtkgl.h>

#ifdef USE_FANCY_CMDLINE
#include <vte/vte.h>
#endif

//#ifdef G_OS_WIN32
//#define WIN32_LEAN_AND_MEAN 1
//#include <windows.h>
//#endif

#include <GL/gl.h>
#include <GL/glu.h>

#include "trackball.h"

#include "qnm.h"

// #defines {{{1
#define MAX_LIGHTS 8
#define MAX_CLIPPLANES 6

#define MAX_SELECTION 300
#define CACHE_LAST 10
#define CACHE_PREDICTED 20
// TODO: heuristic predictions in a separate thread
#define CACHE_SIZE (CACHE_LAST+CACHE_PREDICTED)

#define STARTING_LAMBDA 2.3

#define ANIMATE_THRESHOLD 25.0

#define VIEW_SCALE_MAX 50.0
#define VIEW_SCALE_MIN 0.0625

#define NUM_SHAPES (MAX_SELECTION+5)

#define GLERROR() { \
	GLenum err = glGetError(); \
	while (err != GL_NO_ERROR) { \
		fprintf(stderr, "glError: \"%s\" caught at %s:%d\n", (char *)gluErrorString(err), __FILE__, __LINE__); \
		err = glGetError(); \
	} \
}

#define HARDCODED_FIVE 5
SMART_PARAMETER(UNREAL_FIVE,HARDCODED_FIVE)

#define CMDLINE_LENGTH 0x77

#define DEFAULT_CAMERA_POS_FNAME "/tmp/glshapest.default_camera_pos"


// }}}1

// typedef struct { ... } state; {{{1
typedef struct { /* state */
	// input data
	int s[3], n;			// n = s[0] * s[1] * s[2];
	//bool have_t;
	//bool have_g;
	//bool have_h;
	tree_of_regions *t;
	imatge *g, in[1];
	imatge *vg, vgstatic[3]; // "vec-gradient" (only for visualization)
	collection_of_hypersurfaces *h; // (TODO)

	// visualization state (transformation matrix, pallette, etc)
	float pallette[3][0x100];
	float imgrang;
	//bool show_slice[3]
	bool show_isosurface;
	int projection_type; // 0 = max, 1 = average, 2 = slice
	//int slice_pos[3];
	float isosurface_threshold;
	int use_wireframe;
	bool show_voxel;
	int voxel_pos[3];

	// trackball
	float view_quat_diff[4];
	float view_quat[4];
	float view_scale;
	float trackball_begin_x;
	float trackball_begin_y;
	float trackball_dx;
	float trackball_dy;
	gboolean animate;
	float global_scale;
	float global_displacement[3];
	guint idle_id;

	// current "live" state
	region *current_region; // (inside the TOS e->t!!!)
	surface_triangulation *active_surface;
	surface_triangulation *current_boundary;
	surface_triangulation current_piece[1]; // TODO: pointer to caches
	tr_tree current_ms_tree[1];		// TODO: pointer to caches
	tr_region *current_ms_region;
	surface_triangulation *current_patch;
	bool tree_visibility;
	bool stats_visibility;
	bool dotty_style;
	bool surface_outlines;
	bool tree_updated;
	estat_nfa enfa[1];

	// "mode"
	enum mode {
		MODE_TOS,
		MODE_MSTREE,
		MODE_MALIST,
		MODE_POINT,
		MODE_CMDLINE
	} mode;

	// misc
	float alpha;
	double lambda;
	bool reverse_orientation;
	int mmc_option;

	// tree of regions descriptors
	int number_of_branches;
	int *branches_base;
	int *branches_size;
	int *branches_which;

	// selection
	collection_of_hypersurfaces selection[1];
	bool selection_normals_draw;
	int selection_visibility;
	bool selection_updated;
	bool colored_selection;

	// lights
	bool lights_updated;
	bool lights_fixed_to_camera; // unused yet!
	bool lights_displayed_as_dots;
	bool light_active[MAX_LIGHTS];
	GLfloat light_ambient[MAX_LIGHTS][4];
	GLfloat light_diffuse[MAX_LIGHTS][4];
	GLfloat light_specular[MAX_LIGHTS][4];
	GLfloat light_position[MAX_LIGHTS][4];
	GLfloat light_model_ambient[4];
	GLfloat light_local_view[1];

	// material
	GLfloat material_ambient[4];
	GLfloat material_diffuse[4];
	GLfloat material_specular[4];
	GLfloat material_emission[4];
	GLfloat material_shininess;

	// clip planes (molt lents!, i no sé perquè, cony!)
	bool clipplane_active[MAX_CLIPPLANES];
	GLdouble clipplane_equation[MAX_CLIPPLANES][4];

	// GLdata (display lists, associated with caching)
	GLuint shape_list_base;
	GLuint shape_surface;
	GLuint shape_current;
	GLuint shape_plane_idx;
	GLuint shape_voxeld;

	// caching (TODO)
	int *cached;		// length = t->n;
	int cache_index;	// 0 <= CACHE_SIZE
	// cache of shape boundaries (N last and M predicted)
	surface_triangulation cached_surfaces[CACHE_SIZE];
	// cache of mumford_shah_trees (corresponding to the above)
	tr_tree cached_ms_trees[CACHE_SIZE];
	bool force_recomputation;

	// planes (projection) state
	bool plane_show;
	int plane_current_projection_type;
	imatge planes[3][NUMBER_OF_PROJECTION_TYPES];
	//float plane_projection_position[3];
	bool white_background;
	bool plane_updated;
	bool show_drawing_of_axes;

	// voxel-drawing state (MODE_POINT)
	// "as an example of code organization"
	int voxel_drawing_state;
	imatge voxel_drawing_image[1];
	int voxel_drawing_point[3];

	// cmdline state (MODE_CMDLINE)
	GtkWidget *cmdline_w;
	int cmdline_pos;
	char cmdline_string[CMDLINE_LENGTH];
	enum mode cmdline_oldmode;


	// variable settings
	int unreal_five;	// stride for fast navigation


	// ugly hacks
	GtkWidget *drawing_area; // extremely ugly hack
	GtkWidget *widget;
	GtkWidget *window;
	int hackmode;
	bool strange_colors;
	char *the_hyphen_c_argument;
	// uglier
	tree_of_regions tstatic[1];
	imatge gstatic[1];
	collection_of_hypersurfaces hstatic[1];
	surface_triangulation single_surface[1];
} state;

static bool global_do_not_add_quats = false;

static state *global_state_pointer_hack; // TODO: remove this variable
static state *pick_state(void *data)
{
	state *s = data;
	assert(s);
	if (s != global_state_pointer_hack)
		error("NOT THE STATE %p != %p\n", (void *)s,
				(void *)global_state_pointer_hack);
	return s;
}

static double global_subp_volume = -42;
static double global_subp_area = -42;
static double global_subp_pperim = -42;
static imatge global_metric[1];

// misc state management {{{1
static void start_lights(state *);
static void start_material(state *);
static void start_materialcolors(state *);
static void
init_view (state *e)
{
	e->view_quat[0] = e->view_quat_diff[0] = 0.0;
	e->view_quat[1] = e->view_quat_diff[1] = 0.0;
	e->view_quat[2] = e->view_quat_diff[2] = 0.0;
	e->view_quat[3] = e->view_quat_diff[3] = 1.0;
	e->view_scale = 1.0;

	start_lights(e);
	start_material(e);
	start_materialcolors(e);
	e->lights_fixed_to_camera = false;
	e->lights_displayed_as_dots = false;

	FORI(MAX_CLIPPLANES) e->clipplane_active[i] = false;
	//e->clipplane_active[0] = true;
	e->clipplane_equation[0][0] = 1;
	e->clipplane_equation[0][1] = 1;
	e->clipplane_equation[0][2] = 1;
	e->clipplane_equation[0][3] = 1;
}

static float view_quat_sign(state *e, int axis)
{
	//return -1;
	float *q = e->view_quat, qdis;
	switch(axis) {
	case 0: qdis = q[2]*q[0] + q[1]*q[3]; break;
	case 1: qdis = q[1]*q[2] - q[0]*q[3]; break;
	case 2: qdis = 1 - 2*(q[1]*q[1] + q[0]*q[0]); break;
	default: error("invalid axis %d", axis);
	}
	return qdis;
}

static void init_projection_planes(state *e)
{
	reconstruct(e->in, e->t);
	build_all_statistical_orthogonal_projections(e->planes, e->in);
	float mi, ma;
	get_min_maxf(e->in->grayplane, e->in->mida, &mi, &ma);
	allibera_imatge(e->in);
	e->plane_current_projection_type = 2;
	void normalitza_de(imatge *, float, float);
	FORJ(3) FORI(NUMBER_OF_PROJECTION_TYPES)
	{
		normalitza_de(e->planes[j] + i, mi, ma);
		if (i == 5)
		{
			imatge *p = e->planes[j]+i;
			FORK(p->mida)
				p->grayplane[k] = 1 - p->grayplane[k];
		}
	}
}

static void init_voxel_drawing(state *e)
{
	e->voxel_drawing_state = 0;
	inicialitza_himatge(e->voxel_drawing_image, e->t);
	FORI(3) e->voxel_drawing_point[i] = 0;
}

static void clean_voxel_drawing(state *e)
{
	allibera_imatge(e->voxel_drawing_image);
}

static void init_branches(state *e)
{
	assert(e->t);
	XMALLOC(e->branches_base, e->t->n);
	XMALLOC(e->branches_size, e->t->n);
	XMALLOC(e->branches_which, e->t->n);
	e->number_of_branches = pick_monotone_branches(
			e->branches_base,
		       	e->branches_size,
		       	e->t);
	which_branch_contains_each_region(
			e->branches_which,
			e->branches_base,
			e->branches_size,
			e->number_of_branches,
			e->t);
}

static void clean_branches(state *e)
{
	xfree(e->branches_base);
	xfree(e->branches_size);
	xfree(e->branches_which);
}

SMART_PARAMETER(ENFA_NBINS,256)
SMART_PARAMETER(ENFA_NTESTS,-1)
SMART_PARAMETER(ENFA_EPSILON,1)
SMART_PARAMETER(ENFA_QUANTILE,0.01)

static void init_state(state *e, tree_of_regions *t,
		imatge *g, collection_of_hypersurfaces *h)
{
	global_state_pointer_hack = e;
	e->t = t;
	e->g = g;
	e->h = h;
	assert(e->t);
	assert(e->g);
	//assert(!e->h);
	e->s[0] = e->t->w;
	e->s[1] = e->t->h;
	e->s[2] = e->t->d;
	e->n = e->t->mida;
	if (e->g) assert(e->t->mida == e->g->mida);
	assert(e->t->has_voxels_of_root);
	//t->pinf.x = 10;
	//t->pinf.y = 20;
	//t->pinf.z = 30;
	assert(pertanyP(t, t->root, t->pinf.x, t->pinf.y, t->pinf.z));
	//assert(false);
	//FORI(3) e->show_projection[i] = false;
	e->show_voxel = false;
	e->show_drawing_of_axes = false;
	e->use_wireframe = 0;
	e->cached = NULL;
	init_collection_of_hypersurfaces(e->selection, 3, MAX_SELECTION);
	XMALLOC(e->selection->data.p, 3*sizeof(float)*MAX_SELECTION);
	e->selection_normals_draw = false;
	e->selection->number_of_pieces = 0;
	if (e->h) {
		FORI(e->h->number_of_pieces)
			add_surface_to_collection_deep(e->selection,i+e->h->t);
		free_collection_of_hypersurfaces_as_much_as_possible(e->h);
		e->h = NULL;
	}
	e->selection_visibility = 1;
	e->alpha = 0.5;
	e->reverse_orientation = true;
	e->mmc_option = 0;
	e->force_recomputation = false;
	e->lambda = 10;
	e->tree_visibility = true;
	e->stats_visibility = true;
	e->dotty_style = false;
	e->surface_outlines = true;
	e->plane_show = true;
	e->white_background = false;
	e->colored_selection = false;

	e->mode = MODE_TOS;

	init_view (e);
	e->animate = FALSE;
	e->vg = NULL;

	e->shape_surface = MAX_SELECTION-1;
	e->shape_current = MAX_SELECTION-1;
	e->shape_plane_idx = MAX_SELECTION;
	e->shape_voxeld = MAX_SELECTION-2;
	e->idle_id = 0;
	// cache stuff
	XMALLOC(e->cached, e->t->n); FORI(e->t->n) e->cached[i] = -1;
	FORI(CACHE_SIZE)
	{
		e->cached_surfaces[i].data.i = -1;
		e->cached_surfaces[i].nv = 0;
	}
	e->current_piece->nv = 0;
	e->current_ms_tree->number_of_things = 0;
	e->cache_index = 0;

	e->hackmode = -1;
	e->strange_colors = false;

	distro_histo(&(e->enfa->dac), ENFA_NBINS(), g->grayplane, g->mida);
	//distro_plot(&(e->enfa->dac));
	distro_ac(&(e->enfa->dac));
	//distro_plot(&(e->enfa->dac));
	e->enfa->epsilon = ENFA_EPSILON();
	e->enfa->quantile = ENFA_QUANTILE();
	if (ENFA_NTESTS() < 1)
	{
		//printf("g->mida = %d\n", g->mida);
		//printf("g->mida^1.6666 = %g\n", pow(g->mida,1.6666));
		e->enfa->ntests = 0x100 * pow(g->mida, 1.6666);
		printf("ntests = %lf\n", e->enfa->ntests);
	}
	else
		e->enfa->ntests = ENFA_NTESTS();

	init_projection_planes(e);
	init_voxel_drawing(e);
	init_branches(e);
	e->unreal_five = UNREAL_FIVE();
}

static void clean_state(state *e)
{
	allibera_tor(e->t);
	if (e->g)
		allibera_imatge(e->g);
	FORI(CACHE_SIZE)
	{
		surface_triangulation *s = e->cached_surfaces + i;
		if (s->data.i >= 0)
			free_surface_triangulation(s);
	}
	if (e->current_piece->nv)
		free_surface_triangulation(e->current_piece);
	if (e->current_ms_tree->number_of_things)
		tr_free_tree(e->current_ms_tree);
	xfree(e->cached);
	xfree(e->selection->data.p);
	free_collection_of_hypersurfaces_as_much_as_possible(e->selection);
	if (e->vg) FORI(3) allibera_imatge(e->vg + i); e->vg = NULL;
	distro_free(&(e->enfa->dac));
	FORJ(3)
	FORI(NUMBER_OF_PROJECTION_TYPES) allibera_imatge(e->planes[j] + i);
	clean_voxel_drawing(e);
	clean_branches(e);
}


static void
print_globals(state *e)
{
	printf("\n");
	tree_of_regions *t = e->t;
	region *r = e->current_region;
	printf("region idx,id = %td %d\n", r - t->the_regions, r->id);
	printf("\tvolume = %d, pvol = %d\n", r->volume, r->proper_volume);
	printf("\t\t(subp_volume = %g, subp_area = %g, subp_perim = %g, 1/C = %g, 1/Cg = %g, Cg=%g)\n", global_subp_volume, global_subp_area, global_subp_pperim, global_subp_volume/global_subp_area, global_subp_volume/global_subp_pperim, global_subp_pperim/global_subp_volume);
	if (r->next_sibling)
		printf("\tnext_sibling volume = %d\n", r->next_sibling->volume);
	printf("\n");
}

static void print_active_surface_statistics(state *e)
{
	surface_triangulation *t = e->active_surface;
	assert(t->colors);// if (!t->colors) return;
	statistics_float s;
	statistics_getf(&s, t->colors, t->nv);
	double sigd = significativitat_desolneux(e->enfa, t->nv, s.min);
	//double sigq = significativitat_quantile(e->enfa, t->nv, s.min);
	//printf("\n");
	printf("\tmin      = %g\n", s.min);
	printf("\tmax      = %g\n", s.max);
	printf("\tmedian   = %g\n", s.median);
	printf("\taverage  = %g\n", s.average);
	printf("\tvariance = %g\n", s.variance);
	printf("\tVERTICES = %d\n", t->nv);
	printf("\tsig_des  = %g\n", sigd);
	//printf("\tsig_qua  = %g\n", sigq);
	//printf("\n");
}

static int min_shape_containingB(tree_of_regions *t, imatge *B)
{
	imatge *x = B;
	region *r = NULL;
	assert(B->mida == t->mida);
	FORI(B->mida)
		if (B->grayplane[i] > 0.5)
		{
			region *s = t->smallest_region[i];
			if (!r) r = s;
			r = min_enclosing_region(t, r, s);
		}
	return r - t->the_regions;
}

static void toggle_animation (state *e, GtkWidget *widget);

// }}}1

// OpenGL drawing stuff {{{1

// vector algebra {{{2
static void
crossprodGLGL(GLfloat v1[3], GLfloat v2[3], GLfloat prod[3])
{
	GLfloat p[3];         /* in case prod == v1 or v2 */

	p[0] = v1[1] * v2[2] - v2[1] * v1[2];
	p[1] = v1[2] * v2[0] - v2[2] * v1[0];
	p[2] = v1[0] * v2[1] - v2[0] * v1[1];
	prod[0] = p[0];
	prod[1] = p[1];
	prod[2] = p[2];
}

static void
normalize(GLfloat v[3])
{
	GLfloat d;

	d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	if (d == 0.0) {
		g_warning("normalize: zero length vector");
		v[0] = d = 1.0;
	}
	d = 1 / d;
	v[0] *= d;
	v[1] *= d;
	v[2] *= d;
}

static void
normalize_d(GLdouble v[3])
{
	GLdouble d;

	d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	if (d == 0.0) {
		g_warning("normalize_d: zero length vector");
		v[0] = d = 1.0;
	}
	d = 1 / d;
	v[0] *= d;
	v[1] *= d;
	v[2] *= d;
}


static void compute_vertex_normal(GLfloat *r, state *e, float *v)
{
	float z = e->reverse_orientation ? -1 : 1;
	FORI(3)
		r[i] = z * trilinear_interpolantf(e->vg + i, v[0], v[1], v[2]);
	if (e->hackmode >= 0)
		r[e->hackmode] *= 2;
	normalize(r);
}

#define DIFF3(_a,_b,_c) { \
	(_c)[0] = (_a)[0] - (_b)[0]; \
	(_c)[1] = (_a)[1] - (_b)[1]; \
	(_c)[2] = (_a)[2] - (_b)[2]; \
}

static void compute_triangle_normal(GLfloat *r,
		state *e, surface_triangulation *t, int *tidx)
{
	GLfloat *n[3]; FORI(3) n[i] = t->v[tidx[i]];
	GLfloat q0[3], q1[3];

	DIFF3(n[0], n[1], q0);
	DIFF3(n[1], n[2], q1);
	if (e && e->reverse_orientation)
		crossprodGLGL(q0, q1, q1);
	else
		crossprodGLGL(q1, q0, q1);
	normalize(q1);
	FORI(3) r[i] = q1[i];
}

static void halfpixel_inline(surface_triangulation *t)
{
	return;
	//FORI(t->nv) FORJ(3)
	//	t->v[i][j] -= 0.5;
}

static void compute_surface_normals(float (*o)[3], surface_triangulation *t,
		bool invert)
{
	DXMALLOC(int, tc, t->nv);
	FORI(t->nv)
	{
		tc[i] = 0;
		FORJ(3)
			o[i][j] = 0;
	}
	FORI(t->nt)
       	{
		float n[3];
		compute_triangle_normal(n, NULL, t, t->t[i]);
		FORJ(3)
		{
			int vi = t->t[i][j];
			assert(vi >= 0);
			assert(vi < t->nv);
			tc[vi] += 1;
			FORK(3)
				o[vi][k] += n[k];
		}
	}
	FORI(t->nv)
	{
		//assert(tc[i] > 0);
		if (invert) FORJ(3) o[i][j] *= -1;
		normalize(o[i]);
	}
	xfree(tc);
}

static void surface_normal_motion(surface_triangulation *t, float lambda)
{
	assert(t->n);
	FORI(t->nv)
		FORJ(3)
			t->v[i][j] += lambda * t->n[i][j];
}

// lighting and materials {{{2

static void start_lights(state *e)
{
	// sanity check
	assert(GL_LIGHT1 == GL_LIGHT0 + 1);
	assert(GL_LIGHT2 == GL_LIGHT0 + 2);
	assert(GL_LIGHT3 == GL_LIGHT0 + 3);
	assert(GL_LIGHT4 == GL_LIGHT0 + 4);
	assert(GL_LIGHT5 == GL_LIGHT0 + 5);
	assert(GL_LIGHT6 == GL_LIGHT0 + 6);
	assert(GL_LIGHT7 == GL_LIGHT0 + 7);

	// safe defaults
	FORI(8) e->light_active[i] = false;
	FORI(8) FORJ(4) e->light_ambient[i][j] = j==3?1:0;
	FORI(8) FORJ(4) e->light_diffuse[i][j] = j==3?1:0;
	FORI(8) FORJ(4) e->light_specular[i][j] = 1;
	FORI(8) FORJ(4) e->light_position[i][j] =j==3?1: 0;
	FORJ(4) e->light_model_ambient[j] = j==3?1:0;
	e->light_local_view[0] = 0;

	// fancy values
	GLfloat light_red[] =    { 1,  0,  0,    1};
	GLfloat light_green[] =  { 0,  1,  0,    1};
	GLfloat light_blue[] =   { 0,  0,  1,    1};

	GLfloat light_redx[] =   { 3,  0,  0,    1};
	GLfloat light_greenx[] = { 0,  3,  0,    1};
	GLfloat light_bluex[] =  { 0,  0,  3,    1};

	GLfloat light_none[] =   { 0,  0,  0,    1};
	GLfloat light_full[] =   { 1,  1,  1,    1};

	GLfloat light_po[] =     { 0,  3,  3,    1};
	GLfloat light_mpo[] =    {-3, -3, -3,    1};

	GLfloat light_model_ambient_default[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light_local_view_default[] = {0.0};

	// fill state struct
	FORJ(4) e->light_ambient[0][j] = light_none[j];
	FORJ(4) e->light_diffuse[0][j] = light_red[j];
	FORJ(4) e->light_position[0][j] = light_redx[j];
	e->light_active[0] = true;

	FORJ(4) e->light_ambient[1][j] = light_none[j];
	FORJ(4) e->light_diffuse[1][j] = light_green[j];
	FORJ(4) e->light_position[1][j] = light_greenx[j];
	e->light_active[1] = true;

	FORJ(4) e->light_ambient[2][j] = light_none[j];
	FORJ(4) e->light_diffuse[2][j] = light_blue[j];
	FORJ(4) e->light_position[2][j] = light_bluex[j];
	e->light_active[2] = true;

	FORJ(4) e->light_ambient[3][j] = light_none[j];
	FORJ(4) e->light_diffuse[3][j] = light_full[j];
	FORJ(4) e->light_position[3][j] = light_mpo[j];
	e->light_active[3] = true;

	FORJ(4) e->light_model_ambient[j] = light_model_ambient_default[j];
	e->light_local_view[0] = light_local_view_default[0];
}

static void locate_ligths(state *e)
{
	FORI(8) if (e->light_active[i]) {
		glLightfv(GL_LIGHT0 + i, GL_AMBIENT, e->light_ambient[i]);
		glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, e->light_diffuse[i]);
		glLightfv(GL_LIGHT0 + i, GL_SPECULAR, e->light_specular[i]);
		glLightfv(GL_LIGHT0 + i, GL_POSITION, e->light_position[i]);
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, e->light_model_ambient);
	glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, e->light_local_view);

}

static void display_lights(state *e)
{
	glEnable (GL_LIGHTING);

	FORI(8)
		if (e->light_active[i])
			glEnable(GL_LIGHT0 + i);
}

static void undisplay_lights(state *e)
{
	FORI(8)
		if (e->light_active[i])
			glDisable(GL_LIGHT0 + i);
	glDisable(GL_LIGHTING);
}

static void set_materialc(state *e, GLfloat *col)
{
	glMaterialfv (GL_FRONT, GL_AMBIENT, col);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, col);
	glMaterialfv (GL_FRONT, GL_SPECULAR, col);
	glMaterialf (GL_FRONT, GL_SHININESS, 0.4*128.0);
	glMaterialfv (GL_BACK, GL_AMBIENT, col);
	glMaterialfv (GL_BACK, GL_DIFFUSE, col);
	glMaterialfv (GL_BACK, GL_SPECULAR, col);
	glMaterialf (GL_BACK, GL_SHININESS, 0.4 * 128.0);
}

static void start_material(state *e)
{
	//typedef struct _MaterialProp
	//{
	//	GLfloat ambient[4];
	//	GLfloat diffuse[4];
	//	GLfloat specular[4];
	//	GLfloat shininess;
	//} MaterialProp;
	//MaterialProp mat_silver = {
	//	{0.19225, 0.19225, 0.19225, 1.0},
	//	{0.50754, 0.50754, 0.50754, 1.0},
	//	{0.508273, 0.508273, 0.508273, 1.0},
	//	0.4
	//};
	//MaterialProp mat_specy = {
	//	{0.19225, 0.19225, 0.19225, 1.0},
	//	{0.50754, 0.50754, 0.50754, 1.0},
	//	{0.908273, 0.908273, 0.908273, 1.0},
	//	0.1
	//};
	//MaterialProp *mat_current = &mat_specy;
	GLfloat ambient[] = {0.19225, 0.19225, 0.19225, 1.0};
	GLfloat diffuse[] = {0.50754, 0.50754, 0.50754, 1.0};
	GLfloat specular[] = {0.508273, 0.508273, 0.508273, 1.0};
	GLfloat emission[] = {0, 0, 0, 1.0};
	GLfloat shininess = 0.4;

	FORJ(4) e->material_ambient[j] = ambient[j];
	FORJ(4) e->material_diffuse[j] = diffuse[j];
	FORJ(4) e->material_specular[j] = specular[j];
	FORJ(4) e->material_emission[j] = emission[j];
	e->material_shininess = shininess;
}

static void set_material(state *e)
{
	glMaterialfv (GL_FRONT, GL_AMBIENT, e->material_ambient);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, e->material_diffuse);
	glMaterialfv (GL_FRONT, GL_SPECULAR, e->material_specular);
	glMaterialfv (GL_FRONT, GL_EMISSION, e->material_emission);
	glMaterialf (GL_FRONT, GL_SHININESS, e->material_shininess * 128.0);

	glMaterialfv (GL_BACK, GL_AMBIENT, e->material_ambient);
	glMaterialfv (GL_BACK, GL_DIFFUSE, e->material_diffuse);
	glMaterialfv (GL_BACK, GL_SPECULAR, e->material_specular);
	glMaterialfv (GL_BACK, GL_EMISSION, e->material_emission);
	glMaterialf (GL_BACK, GL_SHININESS, e->material_shininess * 128.0);
}

static void randomize_materialcolors(state *e)
{
	if (!e->colored_selection)
		e->colored_selection = true;
	else {
		float *f = e->selection->data.p;
		assert(f);
		FORI(MAX_SELECTION) FORJ(3) f[3*i+j] = random_uniform_bounds(0.2, 1);
	}
}

static void start_materialcolors(state *e)
{
	float def[6][3] = {
		{.8, .3, .3},
		{.3, .8, .3},
		{.3, .3, .8},
		{.3, .8, .8},
		{.8, .3, .8},
		{.8, .8, .3}
	};
	float *f = e->selection->data.p;
	assert(f);
	FORI(MAX_SELECTION) FORJ(3) f[3*i+j] = random_uniform_bounds(0.2, 1);
	FORI(6) FORJ(3) f[3*i+j] = def[i][j];
}

// dumping primitives {{{2
static void
dump_shaded_triangle(state *e, surface_triangulation *t, int *tidx)
{
	GLfloat n[3];

	glBegin(GL_TRIANGLES);
	if (e->vg)
	{
		FORI(3)
		{
			if (t->n)
				FORK(3) n[k] = t->n[tidx[i]][k];
			else
				compute_vertex_normal(n, e, t->v[tidx[i]]);
			glNormal3fv(n);
			glVertex3fv(t->v[tidx[i]]);
		}
	} else {
		compute_triangle_normal(n, e, t, tidx);
		glNormal3fv(n);
		FORI(3)
			glVertex3fv(t->v[tidx[i]]);
	}
	glEnd();
}

static void
dump_outlined_triangle(state *e, surface_triangulation *t, int *tidx)
{
	glBegin(GL_LINE_LOOP);
	FORI(3)
		glVertex3fv(t->v[tidx[i]]);
	glEnd();
}

static void
dump_colored_triangle(state *e, surface_triangulation *t, int *tidx)
{
	glBegin(GL_TRIANGLES);
	FORI(3)
	{
		float g = t->colors[tidx[i]] / e->imgrang;
		if (e->strange_colors)
		{
			GLfloat d[3], G[3];
			compute_triangle_normal(d, e, t, tidx);
			FORJ(3) {
			float *v = t->v[tidx[j]];
			G[j]=trilinear_interpolantf(e->vg+j,v[0],v[1],v[2]);
			}
			normalize(G);
			g = producte_escalarf(d, G, 3);
			if (e->reverse_orientation) g *= -1;
			g = (g + 1)/2;
			g = 1 - g;
			glColor3f(g>0.5?0.9:0.1, cos(3*g), sin(3*g));
		}
		else
			glColor3f(0.9, cos(3*g), sin(3*g));
		glVertex3fv(t->v[tidx[i]]);
	}
	glEnd();
}

static void draw_pped(GLfloat *a, GLfloat *b, bool solid)
{
	GLenum type = solid ? GL_QUADS : GL_LINE_LOOP;
	GLfloat v[8][3] =
	{
		{a[0], a[1], a[2]}, // 0
		{b[0], a[1], a[2]}, // 1
		{a[0], b[1], a[2]}, // 2
		{b[0], b[1], a[2]}, // 3
		{a[0], a[1], b[2]}, // 4
		{b[0], a[1], b[2]}, // 5
		{a[0], b[1], b[2]}, // 6
		{b[0], b[1], b[2]}  // 7
	};
	GLint f[6][4] =
	{
		{0, 2, 3, 1},
		{1, 3, 7, 5},
		{4, 5, 7, 6},
		{0, 4, 6, 2},
		{0, 1, 5, 4},
		{2, 6, 7, 3}
	};
	GLfloat n[6][3] = 
	{
		{0, 0, -1},
		{1, 0, 0},
		{0, 0, 1},
		{-1, 0, 0},
		{0, -1, 0},
		{0, 1, 0}
	};
	FORI(6)
	{
		glBegin(type);
		glNormal3fv(n[i]);
		FORJ(4)
			glVertex3fv(v[f[i][j]]);
		glEnd();
	}
}

static void draw_pyrq(GLfloat *a, GLfloat *b, GLfloat *c, GLfloat *d,
		GLfloat *v, bool solid)
{
	GLenum t = solid ? GL_TRIANGLES : GL_LINE_LOOP;
	glBegin(t); glVertex3fv(a); glVertex3fv(b); glVertex3fv(v); glEnd();
	glBegin(t); glVertex3fv(b); glVertex3fv(c); glVertex3fv(v); glEnd();
	glBegin(t); glVertex3fv(c); glVertex3fv(d); glVertex3fv(v); glEnd();
	glBegin(t); glVertex3fv(d); glVertex3fv(a); glVertex3fv(v); glEnd();
}

// draw a cool point localizer
#define V(x,y,z) (GLfloat []){x,y,z}
static void draw_sixpyrq(imatge *x, float *pp)
{
	static const float ESLACK = 0.008;
	GLfloat p[3]; FORI(3) p[i] = floor(pp[i]) - 0.5;
	GLfloat q[3]; FORI(3) q[i] = p[i] + 1;
	float W = x->w - 1; float H = x->h - 1; float D = x->d - 1;
	FORL(3) p[l] -= ESLACK;
	FORL(3) q[l] += ESLACK;
	glLineWidth(2.0);
	draw_pped(p, q, false);
	glLineWidth(1.0);
	glColor3f(0,0,1);
	bool so = false;
draw_pyrq(V(p[0],p[1],0),V(q[0],p[1],0),V(q[0],q[1],0),V(p[0],q[1],0),pp,so);
draw_pyrq(V(p[0],p[1],D),V(q[0],p[1],D),V(q[0],q[1],D),V(p[0],q[1],D),pp,so);
draw_pyrq(V(p[0],0,p[2]),V(q[0],0,p[2]),V(q[0],0,q[2]),V(p[0],0,q[2]),pp,so);
draw_pyrq(V(p[0],H,p[2]),V(q[0],H,p[2]),V(q[0],H,q[2]),V(p[0],H,q[2]),pp,so);
draw_pyrq(V(0,p[1],p[2]),V(0,q[1],p[2]),V(0,q[1],q[2]),V(0,p[1],q[2]),pp,so);
draw_pyrq(V(W,p[1],p[2]),V(W,q[1],p[2]),V(W,q[1],q[2]),V(W,p[1],q[2]),pp,so);
}
#undef V

// dumping planes (ugly) {{{2

#define TOLERANCE 0.01
static void dump_x_projection_to_gl(state *e)
{
	imatge *p = e->planes[0] + e->plane_current_projection_type;
	assert(p->d == 1);

	//glPushMatrix();
	//GLfloat *d = e->global_displacement;
	//glTranslatef( d[0], d[1], d[2]);
	//glScalef(e->global_scale, e->global_scale, e->global_scale);
	glBegin(GL_TRIANGLES);
	float D = 0;
	FORJ(p->h)FORI(p->w)
	{
		float g = p->t[0][j][i], h = 0.5;
		assert(g >= 0-TOLERANCE);
		assert(g <= 1+TOLERANCE);
		glColor3f(g, g, g);

		glVertex3f(D, i-h, j-h);
		glVertex3f(D, i-h, j+h);
		glVertex3f(D, i+h, j-h);

		glVertex3f(D, i-h, j+h);
		glVertex3f(D, i+h, j+h);
		glVertex3f(D, i+h, j-h);
	}
	glEnd();
	//glPopMatrix();

	GLERROR();
}

static void dump_y_projection_to_gl(state *e)
{
	imatge *p = e->planes[1] + e->plane_current_projection_type;
	assert(p->d == 1);

	//glPushMatrix();
	//GLfloat *d = e->global_displacement;
	//glTranslatef( d[0], d[1], d[2]);
	//glScalef(e->global_scale, e->global_scale, e->global_scale);
	glBegin(GL_TRIANGLES);
	float D = 0;
	FORJ(p->h)FORI(p->w)
	{
		float g = p->t[0][j][i], h = 0.5;
		assert(g >= 0-TOLERANCE);
		assert(g <= 1+TOLERANCE);
		glColor3f(g, g, g);

		glVertex3f(i-h, D, j-h);
		glVertex3f(i-h, D, j+h);
		glVertex3f(i+h, D, j-h);

		glVertex3f(i-h, D, j+h);
		glVertex3f(i+h, D, j+h);
		glVertex3f(i+h, D, j-h);
	}
	glEnd();
	//glPopMatrix();

	GLERROR();
}

static void dump_z_projection_to_gl(state *e)
{
	imatge *p = e->planes[2] + e->plane_current_projection_type;
	assert(p->d == 1);

	//glPushMatrix();
	//GLfloat *d = e->global_displacement;
	//glTranslatef( d[0], d[1], d[2]);
	//glScalef(e->global_scale, e->global_scale, e->global_scale);
	glBegin(GL_TRIANGLES);
	float D = 0;
	FORJ(p->h)FORI(p->w)
	{
		float g = p->t[0][j][i], h = 0.5;
		assert(g >= 0-TOLERANCE);
		assert(g <= 1+TOLERANCE);
		glColor3f(g, g, g);

		glVertex3f(j-h, i-h, D);
		glVertex3f(j+h, i-h, D);
		glVertex3f(j-h, i+h, D);

		glVertex3f(j+h, i-h, D);
		glVertex3f(j+h, i+h, D);
		glVertex3f(j-h, i+h, D);
	}
	glEnd();
	//glPopMatrix();

	GLERROR();
}
//static void dump_projections_to_gl(state *e)
//{
//	dump_x_projection_to_gl(e);
//	dump_y_projection_to_gl(e);
//}


// dumping structures {{{2

static void
dump_region_dots(state *e)
{
	//glPushMatrix();
	//GLfloat *d = e->global_displacement;
	//glTranslatef( d[0], d[1], d[2]);
	//glScalef(e->global_scale, e->global_scale, e->global_scale);

	region *r = e->current_region;
	assert(r >= e->t->the_regions);
	assert(r < e->t->the_regions + e->t->n);
	voxel_coords *v = r->voxels;
	assert(v);

	glBegin(GL_POINTS);
	glColor3f(1, 0, 0);
	for (int i = 0; i < r->volume; i++, v++)
		glVertex3f(v->x, v->y, v->z);
	glEnd();

	//glPopMatrix();
	GLERROR();
}

// this function is meant to be called to build display lists
static void
dump_surface_triangulation_to_gl(state *e, surface_triangulation *t,
		bool colored, bool wireframed)
{
	assert(sizeof(GLfloat) == sizeof(float));
	colored = colored && t->colors;

	if (!colored)
	{
		if (e->colored_selection) {
		assert(e->selection->t <= t);
		assert(e->selection->number_of_pieces+e->selection->t > t);
		assert(e->selection->data.p);
		int sfidx = t - e->selection->t;
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnable(GL_COLOR_MATERIAL);
		GLfloat *c = e->selection->data.p;
		c += 3*sfidx;
		fprintf(stderr, "sfidx = %d\n", sfidx);
		fprintf(stderr, "c = {%g %g %g}\n", c[0], c[1], c[2]);
		//GLfloat c[3] = {0.3, 0.8, 0.3};
		glColor3fv(c);
		}
		display_lights(e);
		set_material(e);
	} else
		if (e->hackmode != -1)
			fprintf(stderr, "hackmode = %d\n", e->hackmode);


	fprintf(stderr, "TDUMP: %scolored surface of %d vertices (mode %d)\n",
			colored ? "" : "non-", t->nv, (int)e->mode);

	//glPushMatrix();
	//GLfloat *d = e->global_displacement;
	//glTranslatef( d[0], d[1], d[2]);
	//glScalef(e->global_scale, e->global_scale, e->global_scale);

	for (int i = 0; i < t->nt; i++)
		if (colored)
			dump_colored_triangle(e, t, t->t[i]);
		else {
			dump_shaded_triangle(e, t, t->t[i]);
		}

	if (e->surface_outlines) {
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glColor3f(0, 0, 0);
		glLineWidth (3.0);
		for (int i = 0; i < t->nt; i++)
			dump_outlined_triangle(e, t, t->t[i]);
		glPopAttrib();
	}

	//glPopMatrix();

	if (!colored)
	{
		glDisable(GL_COLOR_BUFFER_BIT);
		undisplay_lights(e);
	}

	GLERROR();
}

static void dump_voxeldrawing_list(state *e)
{
	static const float ESLACK = 0.005;

	if (!e->voxel_drawing_state) return;
	glNewList(e->shape_list_base + e->shape_voxeld, GL_COMPILE);

	//glPushMatrix();
	//GLfloat *d = e->global_displacement;
	//glTranslatef( d[0], d[1], d[2]);
	//glScalef(e->global_scale, e->global_scale, e->global_scale);

	if (e->white_background)
		glColor3f(0, 0, 0);
	else
		glColor3f(1, 0, 0.3);
		//glColor3f(1, 0.73, 0);
		//glColor3f(0.8, 0, 0.3);

	{
		GLfloat a[] = {0, 0, 0};
		GLfloat b[] = {e->t->w-1, e->t->h-1, e->t->d-1};
		draw_pped(a, b, false);
	}

	if (e->voxel_drawing_state < 2)
		goto endofu;

	{
	int *p = e->voxel_drawing_point;
	float aa[] = {p[0], p[1], p[2]};
	if (e->voxel_drawing_state > 2)
		draw_sixpyrq(e->voxel_drawing_image, aa);

	imatge *x = e->voxel_drawing_image;
	GLfloat col[] = {1, 0.7, 0.7};
	GLfloat bla[] = {0, 0, 0};

	display_lights(e);
	FORK(x->d)FORJ(x->h)FORI(x->w)
	{
		if (x->t[k][j][i] > 0)
		{
			GLfloat a[] = {i, j, k};
			FORL(3) a[l] -= 0.5;
			GLfloat b[3]; FORL(3) b[l] = a[l] + 1;
			set_materialc(e, col);
			draw_pped(a, b, true);

			FORL(3) a[l] -= ESLACK;
			FORL(3) b[l] += ESLACK;
			set_materialc(e, bla);
			draw_pped(a, b, false);

		}
	}
	undisplay_lights(e);
	}

endofu:
	//glPopMatrix();

	glEndList();
}

// main dumpage {{{2

// dump to lists
static void dump_objects(state *e)
{
	// TODO : Now that I understand the purpose of glLists,
	// clean the whole visibilidy/updatedness mess in this function
	// and the "expose" callback

	e->view_quat_diff[0] = 0.0;
	e->view_quat_diff[1] = 0.0;
	e->view_quat_diff[2] = 0.0;
	e->view_quat_diff[3] = 1.0;
	// dump selection:
	if (e->selection_updated)
	FORI(e->selection->number_of_pieces)
	{
		e->selection_updated = false;
		glNewList(e->shape_list_base + i, GL_COMPILE);
		dump_surface_triangulation_to_gl(e, e->selection->t + i,
				false, false);
		if (e->selection_normals_draw)
		{
			//glPushMatrix();
			//GLfloat *d = e->global_displacement;
			//glTranslatef( d[0], d[1], d[2]);
			//glScalef(e->global_scale, e->global_scale, e->global_scale);

			surface_triangulation *t = e->selection->t + i;
			assert(t->n);
			float nk = 0.5;
			FORJ(t->nv)
			{
				glBegin(GL_LINES);
				float a[3], b[3];
				FORK(3) a[k] = t->v[j][k];
				FORK(3) b[k] = a[k] + nk*t->n[j][k];
				glColor3f(1, 0, 1);
				glVertex3fv(a);
				glVertex3fv(b);
				glEnd();
			}
			//glPopMatrix();
		}
		glEndList();
	}
	// dump planes:
	if (e->plane_show && e->plane_updated)
	{
		e->plane_updated = false;
		glNewList(e->shape_list_base + e->shape_plane_idx, GL_COMPILE);
		dump_x_projection_to_gl(e);
		glEndList();

		glNewList(e->shape_list_base +e->shape_plane_idx+1,GL_COMPILE);
		dump_y_projection_to_gl(e);
		glEndList();

		glNewList(e->shape_list_base +e->shape_plane_idx+2,GL_COMPILE);
		dump_z_projection_to_gl(e);
		glEndList();
	}
	// dump drawing:
	dump_voxeldrawing_list(e);
	// dump shape:
	if (e->tree_visibility && e->tree_updated)
	{
		e->tree_updated = false;
		glNewList(e->shape_list_base + e->shape_surface, GL_COMPILE);
		if (e->dotty_style)
			dump_region_dots(e);
		else
			dump_surface_triangulation_to_gl(e,
					e->active_surface, true, true);
		glEndList();
	}
	//if (!e->voxel_drawing_state)
	//{
	//	//glColor3f(1, 0.67, 0); //orange
	//	if (e->plane_show)
	//		glColor3f(1, 0.73, 0); //orange
	//	else
	//		glColor3f(0, 0, 0);
	//	gdk_gl_draw_cube (FALSE, 2.0);
	//}
}

static void update_global_hacks(surface_triangulation *t)
{
	double enclosed_volume(surface_triangulation*, float[3]);
	double isotropic_riemannian_area(surface_triangulation*,imatge*);
	float p[3] = {0, 0, 0};
	global_subp_volume = enclosed_volume(t, p);
	global_subp_area = isotropic_riemannian_area(t, NULL);
	global_subp_pperim = isotropic_riemannian_area(t, global_metric);
}


// surface cache management {{{1

// pick a surface from the cache if it is there
// otherwise compute it and store at the cache
static surface_triangulation *pick_boundary_triangulation(state *e, region *r)
{
	surface_triangulation *s = NULL;
	int ridx = r - e->t->the_regions;
	assert(ridx >= 0);
	assert(ridx < e->t->n);
	printf("looking for boundary of ridx %d (value=%g) \n", ridx, r->value);
	if (e->cached[ridx] >= 0 && !e->force_recomputation)
	{
		printf("it seems that it is already cached, on position %d\n",
				e->cached[ridx]);
		assert(e->cached[ridx] < CACHE_SIZE);
		s = e->cached_surfaces + e->cached[ridx];
		if (s->data.i != ridx)
			error("%d != %d\n", s->data.i, ridx);
		assert(s->data.i == ridx);
	} else {
		if (!e->force_recomputation)
			assert(-1 == e->cached[ridx]);
		assert(e->cache_index >= 0);
		assert(e->cache_index < CACHE_SIZE);
		s = e->cached_surfaces + e->cache_index;
		if (s->nv)
		{
			free_surface_triangulation(s);
			e->cached[s->data.i] = -1;
		}
		e->cached[ridx] = e->cache_index;
		// TODO: the following cut should be a separate function
		// --begin cut
		find_cleant_triangulation_of_border_alpha(s, e->t, r,
				e->mmc_option, e->alpha);
		fes_connectivity(s);
		orient_consistently(s);
		halfpixel_inline(s);
		if (e->g)
			colorize_vertices(s, e->g);
		// --end cut
		s->data.i = ridx;
		e->cache_index += 1;
		if (e->cache_index == CACHE_SIZE)
			e->cache_index = 0;
	}
	assert(s);
	update_global_hacks(s);
	return s;
}

#if 0
// pick it from the cache if it is there
// otherwise compute it and store at the cache
static tr_tree *pick_surface_mstree(state *e, region *r)
{
	surface_triangulation *s = pick_boundary_triangulation(e, r);
	int ridx = s->data.i;
	assert(ridx == r - e->t->the_regions);
	assert(e->cached[ridx] >= 0);
	assert(e->cached[ridx] < CACHE_SIZE);
	error("not implemented");
}
#endif

#if 0
// meaning of the e->cached array:
// e->cached[i] describes the caching status of the ith region of the tree
// e->cached[i] == -1 means that the ith region is not cached
// e->cached[i] == k, 0<=k<CACHE_SIZE, means the boundary is cached at k
// e->cached[i] == k+CACHE_SIZE, 0<=k<CACHE_SIZE, means the partition tree
// 							is cached at k
// it is not possible to cache the tree withoud the boundary

static int cache_index_surface(state *e, int i)
{
	int k = e->cached[i];
	return (k >= 0 && k < 2*CACHE_SIZE) ? k % CACHE_SIZE : -1;
}

static int cache_index_tree(state *e, int i)
{
	int k = e->cached[i];
	return (k >= CACHE_SIZE && k < 2*CACHE_SIZE) ? k - CACHE_SIZE : -1;
}
#endif

static void predict_cache(state *e)
{
	error("not implemented");
}


// GTK callbacks {{{1

//static void debug_callback(char *s)
//{
//	printf("CBACK: %s called\n");
//}

static gboolean process_key_anywhere(state *, int);

// CALLBACK realize {{{2
static void
realize (GtkWidget *widget, gpointer data)
{
	state *e = pick_state(data);
	printf("CBACK: realize called!\n");
	GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
	GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);


	/*** OpenGL BEGIN ***/
	if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
		return;

	if (e->white_background)
		glClearColor (10, 1, 1, 1);
	else
		glClearColor (0.5, 0.8, 0.5, 1.0);
	glClearDepth (1.0);
	glShadeModel(GL_SMOOTH);


	glFrontFace (GL_CW);

	glEnable (GL_AUTO_NORMAL);
	glEnable (GL_NORMALIZE);
	glEnable (GL_DEPTH_TEST);
	glEnable (GL_POINT_SMOOTH);
	if (false) {
		// fancy lines
		//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glEnable (GL_LINE_SMOOTH);
		glEnable (GL_BLEND);
		//glBlendFunc (GL_SRC_ALPHA, GL_ONE);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glLineWidth (1.0);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

		// fancy polygons
		//glEnable (GL_POLYGON_SMOOTH);
		//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	}
	glDepthFunc (GL_LESS);

	/* Shape display lists */
	e->shape_list_base = glGenLists (NUM_SHAPES);

	e->selection_updated = true;
	e->plane_updated = true;
	e->tree_updated = true;
	dump_objects(e);

	gdk_gl_drawable_gl_end (gldrawable);
	/*** OpenGL END ***/


	return;
}

// to be called once, AFTER, the first expose event
static void process_the_hyphen_c_argument(state *e)
{
	if (e->the_hyphen_c_argument)
	{
		char *s = e->the_hyphen_c_argument;
		printf("treating the -c argument: \"%s\"\n", s);
		while (*s) {
			int k = -42;
			if (*s == '\\') {
				switch (s[1]) {
				case 'n': k = '\n'; s += 1; break;
				default:  k = s[1]; s += 1; break;
				}
			} else
				k = *s;
			assert(k >= 0);
			process_key_anywhere(e, k);
			s += 1;
		}
		e->the_hyphen_c_argument = NULL;
	}//else error("do not process the hypen c argument twice!");
}

// CALLBACK configure_event {{{2
static gboolean
configure_event (GtkWidget         *widget,
		GdkEventConfigure *event,
		gpointer           data)
{
	printf("CBACK: configure (reshape) called!\n");
	GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
	GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

	GLfloat w = widget->allocation.width;
	GLfloat h = widget->allocation.height;
	GLfloat aspect;

	/*** OpenGL BEGIN ***/
	if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
		return FALSE;

	glViewport (0, 0, w, h);

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	if (w > h)
	{
		aspect = w / h;
		glFrustum (-aspect, aspect, -1.0, 1.0, 5.0, 60.0);
	}
	else
	{
		aspect = h / w;
		glFrustum (-1.0, 1.0, -aspect, aspect, 5.0, 60.0);
	}

	glMatrixMode (GL_MODELVIEW);

	gdk_gl_drawable_gl_end (gldrawable);
	/*** OpenGL END ***/

	return TRUE;
}


// CALLBACK expose_event {{{2
static gboolean
expose_event (GtkWidget      *widget,
		GdkEventExpose *event,
		gpointer        data)
{
	state *e = pick_state(data);
	widget = e->drawing_area;
	static int expocount = 0;
	//printf("CBACK: expose called! (%d)\n", expocount++);
	GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
	GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

	float m[4][4];

	/*** OpenGL BEGIN ***/
	if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
		return FALSE;

	if (e->white_background)
		glClearColor (1, 1, 1, 1);
	else
		glClearColor (0.5, 0.8, 0.5, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity ();

	/* View transformation. */
	glTranslatef (0.0, 0.0, -10.0);


	glScalef (e->view_scale, e->view_scale, e->view_scale);
	if (!global_do_not_add_quats)
		add_quats (e->view_quat_diff, e->view_quat, e->view_quat);
	build_rotmatrix (m, e->view_quat);
	glMultMatrixf (&m[0][0]);

	/* set ligths (but not lighting!) */
	locate_ligths(e);

	if (true && e->clipplane_active[0]) {
		GLdouble *ce = e->clipplane_equation[0];
		GLdouble x[3] = {0, 0, 0};
		GLdouble y[3]; FORJ(3) y[j] = ce[j];
		normalize_d(y);
		GLdouble q[4]; FORJ(3) q[j] = -y[j];
		FORJ(3) y[j] *= ce[3];
		q[3] = ce[3];

		glColor3f(0, 0, 0);

		glPointSize(5);
		glBegin(GL_POINTS);
		glVertex3dv(x);
		glEnd();

		if (true) {
		glBegin(GL_LINES);
		glVertex3dv(x);
		glVertex3dv(y);
		glEnd();
		}

		glEnable(GL_CLIP_PLANE0);
		glClipPlane(GL_CLIP_PLANE0, q);

	}

	glPushMatrix();
	GLfloat *d = e->global_displacement;
	glTranslatef( d[0], d[1], d[2]);
	glScalef(e->global_scale, e->global_scale, e->global_scale);


	/* Render shapes */
	if (e->selection_visibility)
		FORI(e->selection->number_of_pieces)
			glCallList(e->shape_list_base + i);
	if (e->voxel_drawing_state)
		glCallList(e->shape_list_base + e->shape_voxeld);
	if (e->tree_visibility)
		glCallList(e->shape_list_base + e->shape_current);
	if (true && e->plane_show)
	{
		static const float ESLACK = 0.04;
		float t;

		glPushMatrix();
		t = view_quat_sign(e, 0) > 0 ? -ESLACK : e->t->w-1 + ESLACK;
		if (e->plane_current_projection_type == 5) t = 1;
		glTranslatef(t, 0, 0);
		glCallList(e->shape_list_base + e->shape_plane_idx);
		glPopMatrix();

		glPushMatrix();
		t = view_quat_sign(e, 1) > 0 ? -ESLACK : e->t->h-1 + ESLACK;
		if (e->plane_current_projection_type == 5) t = 1;
		glTranslatef(0, t, 0);
		glCallList(e->shape_list_base + e->shape_plane_idx+1);
		glPopMatrix();

		glPushMatrix();
		t = view_quat_sign(e, 2) > 0 ? -ESLACK : e->t->d-1 + ESLACK;
		if (e->plane_current_projection_type == 5) t = 1;
		glTranslatef(0, 0, t);
		glCallList(e->shape_list_base + e->shape_plane_idx+2);
		glPopMatrix();
	}
	glPopMatrix();
	glDisable(GL_CLIP_PLANE0);
	if (e->lights_displayed_as_dots)
	{
		glPointSize(13);
		FORI(MAX_LIGHTS)
			if (e->light_active[i])
			{
				float *x = e->light_position[i];
				if (x[3] == 1) {
					glBegin(GL_POINTS);
					glColor4fv(e->light_diffuse[i]);
					glVertex3fv(x);
					glEnd();
				} else if (x[3] == 0) {
					float y[3]; FORJ(3) y[j] = 100*x[j];
					glBegin(GL_LINES);
					glColor4fv(e->light_diffuse[i]);
					glVertex3fv(x);
					glVertex3fv(y);
					glEnd();
				} else error("BAD w=%g", x[3]);
			}
	}
	if (e->show_drawing_of_axes) {
		GLfloat x[4][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
		FORI(3) {
			glBegin(GL_LINES);
			glColor3fv(x[i]);
			glVertex3fv(x[3]);
			glVertex3fv(x[i]);
			glEnd();
		}
	}



	GLERROR();

	/* Swap buffers */
	if (gdk_gl_drawable_is_double_buffered (gldrawable))
		gdk_gl_drawable_swap_buffers (gldrawable);
	else
		glFlush ();

	gdk_gl_drawable_gl_end (gldrawable);
	/*** OpenGL END ***/

	static bool ran_hyphen_c = false;
	if (!ran_hyphen_c) {
		ran_hyphen_c = true;
		process_the_hyphen_c_argument(e);
	}


	return TRUE;
}


// CALLBACK button_press_event {{{2
static gboolean
button_press_event (GtkWidget      *widget,
		GdkEventButton *event,
		gpointer        data)
{
	state *e = pick_state(data);
	printf("CBACK: button_press called! (%g %g)[%u %u] (%g %g)\n",
			event->x, event->y, event->state, event->button,
			event->x_root, event->y_root);
	if (e->animate)
	{
		if (event->button == 1)
			toggle_animation (e, widget);
	}
	else
	{
		e->view_quat_diff[0] = 0.0;
		e->view_quat_diff[1] = 0.0;
		e->view_quat_diff[2] = 0.0;
		e->view_quat_diff[3] = 1.0;
	}

	e->trackball_begin_x = event->x;
	e->trackball_begin_y = event->y;

	return FALSE;
}

// CALLBACK button_release_event {{{2
static gboolean
button_release_event (GtkWidget      *widget,
		GdkEventButton *event,
		gpointer        data)
{
	state *e = pick_state(data);
	printf("CBACK: button_release called! [%u]\n", event->button);
	float dx = e->trackball_dx;
	float dy = e->trackball_dy;
	if (!e->animate)
	{
		if (event->button == 1 &&
				((dx*dx + dy*dy) > ANIMATE_THRESHOLD))
			toggle_animation (e, widget);
	}

	e->trackball_dx = 0.0;
	e->trackball_dy = 0.0;

	return FALSE;
}

// CALLBACK motion_notify_event {{{2
static gboolean
motion_notify_event (GtkWidget      *widget,
		     GdkEventMotion *event,
		     gpointer        data)
{
	state *e = pick_state(data);
	//printf("CBACK: motion_notify called! (%g %g) %u (%g %g)\n",
	//		event->x, event->y, event->state,
	//		event->x_root, event->y_root);
	float w = widget->allocation.width;
	float h = widget->allocation.height;
	float x = event->x;
	float y = event->y;
	gboolean redraw = FALSE;

	/* Rotation. */
	if (event->state & GDK_BUTTON1_MASK)
	{
		trackball (e->view_quat_diff,
				(2.0 * e->trackball_begin_x - w) / w,
				(h - 2.0 * e->trackball_begin_y) / h,
				(2.0 * x - w) / w,
				(h - 2.0 * y) / h);

		e->trackball_dx = x - e->trackball_begin_x;
		e->trackball_dy = y - e->trackball_begin_y;

		redraw = TRUE;
	}

	/* Scaling. */
	if (event->state & GDK_BUTTON2_MASK)
	{
		e->view_scale = e->view_scale *
					(1.0 + (y - e->trackball_begin_y) / h);
		if (e->view_scale > VIEW_SCALE_MAX)
			e->view_scale = VIEW_SCALE_MAX;
		else if (e->view_scale < VIEW_SCALE_MIN)
			e->view_scale = VIEW_SCALE_MIN;

		redraw = TRUE;
	}

	e->trackball_begin_x = x;
	e->trackball_begin_y = y;

	if (redraw && !e->animate)
		gdk_window_invalidate_rect (widget->window,
				&widget->allocation, FALSE);

	return TRUE;
}

static int pick_corresponding_char(int k)
{
	int r = -42;
	if (k>= 0 && k < 0x100)	r = k;
	//else if (k == GDK_Escape)	r = 't'; // vi-like
	//else if (k == GDK_Escape)	r = 'q'; // (nenazas)
	else if (k == GDK_Escape)	r = '\e';
	else if (k == GDK_Return)	r = '\n';
	else if (k == GDK_Up)		r = 'k';
	else if (k == GDK_Down)		r = 'j';
	else if (k == GDK_Left)		r = 'h';
	else if (k == GDK_Right)	r = 'l';
	else if (k == GDK_BackSpace)	r = '\b';
	else if (k == GDK_Page_Up)	r = '\v';
	else if (k == GDK_Page_Down)	r = '\f';
	else r = 0;
	assert(r >= 0 && r < 0x100);
	return r;
}

// CALLBACK key_press_event {{{2
static void cmdline_next_key(state *, int);
static void (* const global_static_table_of_key_functions[][0x100])(state *);

static gboolean process_key_anywhere(state *e, int key)
{
	if (e->mode == MODE_CMDLINE)
		cmdline_next_key(e, key);

	void (*const*f)(state*) = global_static_table_of_key_functions[e->mode];
	f += key;
	if (*f)
	{
		(*f)(e);
		return TRUE;
	} else return FALSE;
}

static gboolean
key_press_event (GtkWidget   *widget,
		 GdkEventKey *event,
		 gpointer     data)
{
	assert(event);
	state *e = pick_state(data);
	e->widget = widget;
	//printf("CBACK: key_press called (%p, %u{%x},'%c')!\n", (void *)e, event->keyval, event->keyval, (char)event->keyval);

	int key = pick_corresponding_char(event->keyval);
	return process_key_anywhere(e, key);
}

// idling {{{2
static gboolean
idle (GtkWidget *widget)
{
	//printf("idle called!\n");
	/* Invalidate the whole window. */
	gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);

	/* Update synchronously. */
	gdk_window_process_updates (widget->window, FALSE);

	return TRUE;
}

//static guint idle_id = 0;

static void
idle_add (state *e, GtkWidget *widget)
{
	printf("idle_add called!\n");
	if (e->idle_id == 0)
	{
		e->idle_id = g_idle_add_full (GDK_PRIORITY_REDRAW,
				(GSourceFunc) idle,
				widget,
				NULL);
	}
}

static void
idle_remove (state *e, GtkWidget *widget)
{
	printf("idle_remove called!\n");
	if (e->idle_id != 0)
	{
		g_source_remove (e->idle_id);
		e->idle_id = 0;
	}
}

// CALLBACK map_event {{{2
static gboolean
map_event (GtkWidget *widget,
	   GdkEvent  *event,
	   gpointer   data)
{
	state *e = pick_state(data);
	printf("CBACK: map_event called!\n");
	if (e->animate)
		idle_add (e, widget);

	return TRUE;
}

// CALLBACK unmap_event {{{2
static gboolean
unmap_event (GtkWidget *widget,
	     GdkEvent  *event,
	     gpointer   data)
{
	state *e = pick_state(data);
	printf("CBACK: unmap_event called!\n");
	idle_remove (e, widget);

	return TRUE;
}

// CALLBACK visibility_notify_event {{{2
static gboolean
visibility_notify_event (GtkWidget          *widget,
			 GdkEventVisibility *event,
			 gpointer            data)
{
	state *e = pick_state(data);
	printf("CBACK: visibility_notify called!\n");
	if (e->animate)
	{
		if (event->state == GDK_VISIBILITY_FULLY_OBSCURED)
			idle_remove (e, widget);
		else
			idle_add (e, widget);
	}

	return TRUE;
}


// synthetic callbacks {{{1

static void
toggle_animation (state *e, GtkWidget *widget)
{
	e->animate = !e->animate;

	if (e->animate)
	{
		idle_add (e, widget);
	}
	else
	{
		idle_remove (e, widget);
		e->view_quat_diff[0] = 0.0;
		e->view_quat_diff[1] = 0.0;
		e->view_quat_diff[2] = 0.0;
		e->view_quat_diff[3] = 1.0;
		gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);
	}
}

static void
treat_new_tos_region(state *e)
{
	/*assert(e && e->mode == MODE_TOS);*/
	//free_surface_triangulation(e->current_surface);
	//produce_triangulation();
	if (!e->dotty_style)
	{
		e->current_boundary = pick_boundary_triangulation(e, e->current_region);
		e->active_surface = e->current_boundary;
	}
	e->tree_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
	print_globals(e);
	if (e->stats_visibility)
		print_active_surface_statistics(e);
}

static void
treat_new_mstree_region(state *e)
{
	assert(e && (e->mode == MODE_MSTREE || (e->mode==MODE_CMDLINE && e->cmdline_oldmode==MODE_MSTREE)));
	fprintf(stderr, "treating new mstree region %td\n",
		e->current_ms_region - e->current_ms_tree->the_regions);
	if (e->current_piece->nv)
		free_surface_triangulation(e->current_piece);
	assert(e->current_piece->nv == 0);
	ms_pick_surface_piece(e->current_piece, e->current_boundary,
			e->current_ms_tree, e->current_ms_region);
	e->active_surface = e->current_piece;
	e->tree_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
	if (e->stats_visibility)
		print_active_surface_statistics(e);
}

static void compute_selection_normals(state *e)
{
	//FORI(MAX_SELECTION) assert(!e->selection->t[i].n);
	FORI(e->selection->number_of_pieces)
	{
		surface_triangulation *t = e->selection->t + i;
		XMALLOC(t->n, t->nv);
		//FORJ(t->nv) FORK(3) e->selection_normals[i][j][k] = 0.5;
		compute_surface_normals(t->n, t, e->reverse_orientation);
		fprintf(stderr, "computed %dth surface normals\n", i);
	}
	fprintf(stderr, "computed selection normals\n");
}

// }}}1

// paremeterless (key) commands {{{1

/********************************
 * START TABLE OF KEY FUNCTIONS *
 ********************************/

static void kf_quit(state *e) // {{{2
{
	gtk_main_quit ();
}

static void kf_reverse_orientation(state *e) // {{{2
{
	e->reverse_orientation = !e->reverse_orientation;
	e->selection_updated = true;
	if (e->selection->t->n)
		FORI(e->selection->number_of_pieces)
			FORJ(e->selection->t[i].nv)
				FORK(3)
					e->selection->t[i].n[j][k] *= -1;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_wireframe(state *e) // {{{2
{
	e->use_wireframe = 2;
	e->tree_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_white_background(state *e) // {{{2
{
	e->white_background = !e->white_background;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}


static void kf_tos_root(state *e) // {{{2
{
	e->current_region = e->t->root->first_child;
	treat_new_tos_region(e);
}

static void kf_mstree_root(state *e) // {{{2
{
	e->current_ms_region = e->current_ms_tree->root;
	treat_new_mstree_region(e);
}

static void kf_tos_daughter(state *e) // {{{2
{
	region *r = e->current_region;
	if (r && r->first_child)
	{
		e->current_region = r->first_child;
		//assert(r != r->first_child); // OH FUCK!
		treat_new_tos_region(e);
	} else {
		printf("this region has no child!\n");
	}
}

static void kf_mstree_daughter(state *e) // {{{2
{
	tr_region *r = e->current_ms_region;
	if (r && r->daughter)
	{
		e->current_ms_region = r->daughter;
		//assert(r != r->first_child); // OH FUCK!
		treat_new_mstree_region(e);
	} else {
		printf("this region has no child!\n");
	}
}

static void kf_tos_daughter_fast(state *e) // {{{2
{
	int cx = e->unreal_five;
	region *r = e->current_region;
	while (cx--)
		if (r && r->first_child)
			r = r->first_child;
	e->current_region = r;
	treat_new_tos_region(e);
}

static void kf_tos_bottom_of_branch(state *e) // {{{2
{
	region *r = e->current_region;
	if (r && r->first_child)
	{
		do
			r = r->first_child;
		while (1 == r->number_of_children);
		e->current_region = r;
		treat_new_tos_region(e);
	} else {
		printf("this region has no child!\n");
	}
}

static void kf_tos_top_of_branch(state *e) // {{{2
{
	printf("WARNING: top of branch loosely implemented (minus-one!)");
	region *r = e->current_region;
	if (r && r->parent)
	{
		while (r->parent->first_child == r && !r->next_sibling)
			r = r->parent;
		e->current_region = r;
		treat_new_tos_region(e);
	} else {
		printf("this region has no parent!\n");
	}
}

static void kf_tos_mother(state *e) // {{{2
{
	region *r = e->current_region;
	if (r->parent != e->t->root)
	{
		e->current_region = r->parent;
		treat_new_tos_region(e);
	} else {
		printf("the parent of this region is root!\n");
	}
}

static void kf_mstree_mother(state *e) // {{{2
{
	tr_region *r = e->current_ms_region;
	if (r->mother != e->current_ms_tree->root)
	{
		e->current_ms_region = r->mother;
		treat_new_mstree_region(e);
	} else {
		printf("the parent of this region is root!\n");
	}
}

static void kf_tos_mother_fast(state *e) // {{{2
{
	int cx = e->unreal_five;
	region *r = e->current_region;
	while (cx--)
		if (r->parent != e->t->root)
			r = r->parent;
	e->current_region = r;
	treat_new_tos_region(e);
}

static void kf_tos_sister(state *e) // {{{2
{
	region *r = e->current_region;
	if (r->next_sibling)
	{
		e->current_region = r->next_sibling;
		treat_new_tos_region(e);
	} else {
		printf("there are no more siblings!\n");
	}
}

static void kf_mstree_sister(state *e) // {{{2
{
	tr_region *r = e->current_ms_region;
	if (r->sister)
	{
		e->current_ms_region = r->sister;
		treat_new_mstree_region(e);
	} else {
		printf("there are no more siblings!\n");
	}
}

static void kf_tos_sister_back(state *e) // {{{2
{
	region *r = e->current_region;
	if (r->parent->first_child == r)
	{
		printf("this is already the first sibling!\n");
	} else {
		region *s = r->parent->first_child;
		while (s->next_sibling != r)
			s = s->next_sibling;
		e->current_region = s;
		treat_new_tos_region(e);
	}
}


static void kf_use_wireframe(state *e) // {{{2
{
	if (0 == e->use_wireframe) e->use_wireframe = 1;
	else if (1 == e->use_wireframe) e->use_wireframe = 0;
	else if (2 == e->use_wireframe) e->use_wireframe = 0;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_selection_visibility(state *e) // {{{2
{
	e->selection_visibility = !e->selection_visibility;
	e->selection_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_tree_visibility(state *e) // {{{2
{
	e->tree_visibility = !e->tree_visibility;
	e->tree_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_dotty_style(state *e) // {{{2
{
	e->dotty_style = !e->dotty_style;
	if (!e->dotty_style) e->force_recomputation = true;
	e->tree_updated = true;
	dump_objects(e);
	if (!e->dotty_style) e->force_recomputation = false;
	expose_event(e->widget, NULL, e);
}

static void kf_increase_gamma(state *e) // {{{2
{
	e->imgrang *= 1.2;
	e->tree_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_decrease_gamma(state *e) // {{{2
{
	e->imgrang /= 1.2;
	e->tree_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_yank_active_surface(state *e) // {{{2
{
	assert(e->active_surface);
	printf("yanking surface %p\n", e->active_surface);
	add_surface_to_collection_deep(e->selection, e->active_surface);
	e->selection_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_delete_last_selection(state *e) // {{{2
{
	if (e->selection->number_of_pieces > 0)
	{
		remove_last_surface_from_collection_deep(e->selection);
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	}
}

static void kf_save_selection_coh(state *e) // {{{2
{
	static const char * const fname = "glshapest.out.coh";
	printf("saving collection of surfaces to file \"%s\"\n", fname);
	save_collection_of_hypersurfaces(e->selection, fname);
}

static void kf_tos_save_tc(state *e) // {{{2
{
	static const char * const fname = "glshapest.out.tc";
	printf("saving tree of shapes to file \"%s\"\n", fname);
	tree_of_regions t[1];
	cprm_tor(t, e->t);
	desa_tor_complet(t, fname);
	allibera_tor(t);
}

static void kf_tos_delete_current_shape(state *e) // {{{2
{
	setflag(e->current_region, REGION_IS_REMOVED);
}

static void markdownrec(region *r, region_flag f)
{
	setflag(r, f);
	region *s = r->first_child;
	while (s)
	{
		markdownrec(s, f);
		s = s->next_sibling;
	}
}
static void kf_tos_delete_current_shape_downwards(state *e) // {{{2
{
	markdownrec(e->current_region, REGION_IS_REMOVED);
}

static void kf_mode_tos(state *e) // {{{2
{
	e->mode = MODE_TOS;
	treat_new_tos_region(e);
}

// static void kf_mode_mstree(state *e) // {{{2
SMART_PARAMETER(SLAM,10) // starting lambda
SMART_PARAMETER_INT(USE_MSULT,0)

static void kf_mode_mstree(state *e)
{
	e->mode = MODE_MSTREE;
	surface_triangulation *s
		= pick_boundary_triangulation(e, e->current_region);
	assert(s->colors);
	ms_partition p[1]; p->data = NULL;
	ms_build_partition_from_triangulated_surface(p, s);
	tr_tree *t = e->current_ms_tree;
	if (USE_MSULT())
		trms_build_ult(t, p);
	else
		ms_build_tree_of_mergings(NULL, t, p, SLAM(), _QNM_INFINITY);
	ms_free_partition_structure(p);
	tr_reorder_children_by_volume_recursively(t, t->root);
	assert(t->root);
	e->current_ms_region = t->root->daughter;
	if (!e->current_ms_region) e->current_ms_region = t->root;
	treat_new_mstree_region(e);
}

static void kf_mode_malist(state *e) { e->mode = MODE_MALIST; }

//static void kf_mstree_root(state *e) {}
//static void kf_mstree_daughter(state *e) {}
//static void kf_mstree_mother(state *e) {}
//static void kf_mstree_sister(state *e) {}
static void kf_mstree_sister_back(state *e) {}

static void kf_compute_nice_gradient(state *e) // {{{2
{
	if (e->vg)
	{
		fprintf(stderr, "setting flat shading\n");
		assert(e->vg == e->vgstatic);
		FORI(3)
			allibera_imatge(e->vg + i);
		e->vg = NULL;
	} else {
		fprintf(stderr, "setting Gouraud shading\n");
		e->vg = e->vgstatic;
		imatge x[1]; reconstruct(x, e->t);
		fes_imatges_de_gradients_centrats(e->vg, e->vg+1, e->vg+2, x,
				perllonga_constant_lleig);
		allibera_imatge(x);
	}
	e->selection_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_surface_outlines(state *e) // {{{2
{
	e->surface_outlines = !e->surface_outlines;
	e->selection_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_colorize_strangely(state *e) // {{{2
{
	if (e->strange_colors)
		e->strange_colors = false;
	else {
		if (!e->vg)
			kf_compute_nice_gradient(e);
		e->strange_colors = true;
	}
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_tos_option_zero(state *e) // {{{2
{
	e->mmc_option = 0;
	fprintf(stderr, "new mmc option = %d\n", e->mmc_option);
	e->force_recomputation = true;
	treat_new_tos_region(e);
	e->force_recomputation = false;
}

static void kf_tos_option_one(state *e) // {{{2
{
	e->mmc_option = 1;
	fprintf(stderr, "new mmc option = %d\n", e->mmc_option);
	e->force_recomputation = true;
	treat_new_tos_region(e);
	e->force_recomputation = false;
}

static void kf_tos_option_two(state *e) // {{{2
{
	e->mmc_option = 2;
	fprintf(stderr, "new mmc option = %d\n", e->mmc_option);
	e->force_recomputation = true;
	treat_new_tos_region(e);
	e->force_recomputation = false;
}

static void kf_rotate_projection_type(state *e) // {{{2
{
	e->plane_current_projection_type += 1;
	e->plane_current_projection_type %= NUMBER_OF_PROJECTION_TYPES;
	fprintf(stderr, "ptype = %d (", e->plane_current_projection_type);
	switch(e->plane_current_projection_type) {
	case 0: fprintf(stderr, "min)\n"); break;
	case 1: fprintf(stderr, "max)\n"); break;
	case 2: fprintf(stderr, "avg)\n"); break;
	case 3: fprintf(stderr, "med)\n"); break;
	case 4: fprintf(stderr, "sam)\n"); break;
	case 5: fprintf(stderr, "prx)\n"); break;
	default: error("ivp");
	}
	e->plane_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_randomize_materialcolors(state *e) // {{{2
{
	randomize_materialcolors(e);
	e->selection_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_toggle_planes(state *e) // {{{2
{
	e->plane_show = !e->plane_show;
	e->plane_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_rotate_hackmode(state *e) // {{{2
{
	switch(e->hackmode)
	{
	case 0: e->hackmode = 1; break;
	case 1: e->hackmode = 2; break;
	case 2: e->hackmode = -1; break;
	case -1: e->hackmode = 0; break;
	default: error("invalid hackmode %d", e->hackmode);
	}
	e->selection_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_print_stats(state *e) // {{{2
{
	print_active_surface_statistics(e);
}

static void kf_toggle_stats(state *e) // {{{2
{
	e->stats_visibility = !e->stats_visibility;
	print_active_surface_statistics(e);
}

static void kf_rotate_point_visibility(state *e) // {{{2
{
	e->voxel_drawing_state += 1;
	e->voxel_drawing_state %= 5;
	switch(e->voxel_drawing_state)
	{
	case 0: fprintf(stderr, "point visibility: INVISIBLE\n"); break;
	case 1: fprintf(stderr, "point visibility: MARKED CUBES\n"); break;
	case 2: fprintf(stderr, "point visibility: UNMARKED CUBES\n"); break;
	case 3: fprintf(stderr, "point visibility: MARKED DOTS\n"); break;
	case 4: fprintf(stderr, "point visibility: UNMARKED DOTS\n"); break;
	default: assert(false);
	}
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_mode_point(state *e) // {{{2
{
	e->mode = MODE_POINT;
	e->voxel_drawing_state = 1;
	kf_rotate_point_visibility(e);
}


// static void kf_point_dir_DIR(state *e) // {{{2
static void pno(imatge *x, int *p)
{
	FORI(3) if(p[i] < 0) p[i] = 0;
	if (p[0] >= x->w) p[0] = x->w - 1;
	if (p[1] >= x->w) p[1] = x->h - 1;
	if (p[2] >= x->w) p[2] = x->d - 1;
}
#define PNOX pno(e->voxel_drawing_image, e->voxel_drawing_point);\
	dump_objects(e); expose_event(e->widget, NULL, e)
static void kf_point_dir_xm(state *e) { e->voxel_drawing_point[0]--; PNOX; }
static void kf_point_dir_xp(state *e) { e->voxel_drawing_point[0]++; PNOX; }
static void kf_point_dir_ym(state *e) { e->voxel_drawing_point[1]--; PNOX; }
static void kf_point_dir_yp(state *e) { e->voxel_drawing_point[1]++; PNOX; }
static void kf_point_dir_zm(state *e) { e->voxel_drawing_point[2]--; PNOX; }
static void kf_point_dir_zp(state *e) { e->voxel_drawing_point[2]++; PNOX; }
#undef PNOX

static void kf_point_put(state *e) // {{{2
{
	int *p = e->voxel_drawing_point;
	e->voxel_drawing_image->t[p[2]][p[1]][p[0]] += 1;
	dump_objects(e);
	expose_event(e->widget, NULL, e);

}

static void kf_point_unput(state *e) // {{{2
{
	int *p = e->voxel_drawing_point;
	e->voxel_drawing_image->t[p[2]][p[1]][p[0]] -= 1;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_point_make_center(state *e) // {{{2
{
	e->global_displacement[0] = 0;//- (e->t->w - 1);
	e->global_displacement[1] = 0;//- (e->t->h - 1);
	e->global_displacement[2] = 0;//- (e->t->d - 1);
	FORI(3) e->global_displacement[i] -= e->voxel_drawing_point[i];
	FORI(3) e->global_displacement[i] *= e->global_scale;
	e->tree_updated = true;
	e->selection_updated = true;
	e->plane_updated = true;
	dump_objects(e);
	expose_event(e->widget, NULL, e);
}

static void kf_point_save_asc(state *e) // {{{2
{
	static const char * const fname = "glshapest.out.asc";
	printf("saving edited image to file \"%s\"\n", fname);
	desa_imatge_general(e->voxel_drawing_image, fname);
}

static void cmdline_next_key(state *, int);
static void kf_mode_cmdline(state *e) // {{{2
{
	gtk_widget_show(e->cmdline_w);
	e->cmdline_oldmode = e->mode;
	e->mode = MODE_CMDLINE;
	e->cmdline_string[0] = 0;
	e->cmdline_pos = 0;
	cmdline_next_key(e, ':');
}

// commands with parameters {{{1

static void save_screenshot(state *e, char *fname)
{
	int w = e->drawing_area->allocation.width;
	int h = e->drawing_area->allocation.height;
	imatge_rgba x[2], *y = x+1;
	inicialitza_imatge_rgba(x, w, h, 1);
	glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, x->grayplane);
	GLERROR();
	inicialitza_imatge_rgba(y, w, h, 1);
	FORI(w)FORJ(h) y->t[0][j][i] = x->t[0][h-j-1][i];
	desa_imatge_ppm(y, fname);
	FORI(2) allibera_imatge_rgba(x+i);
}

static void save_camera_position_f(state *e, FILE *f)
{
	FORI(4) fprintf(f, "%g ", e->view_quat[i]);
	fprintf(f, "%g\n", e->view_scale);
}
static void save_camera_position(state *e, char *fname)
{ FILE *f = xfopen(fname, "w"); save_camera_position_f(e, f); xfclose(f); }

static bool set_camera_position(state *e, char *posstring)
{
	float t[5] = {-42, -42, -42, -42, -42};
	fprintf(stderr, "SETTING CAMERA POSITION FROM \"%s\"\n", posstring);
	int r = sscanf(posstring, "%g %g %g %g %g", t, t+1, t+2, t+3, t+4);
	if (r == 5) {
		fprintf(stderr, "set camera position (%g %g %g %g; %g)\n",
				t[0], t[1], t[2], t[3], t[4]);
		FORI(4) e->view_quat[i] = t[i];
		e->view_scale = t[4];
		//expose_event(e->widget, NULL, e);
		return true;
	} else {
		fprintf(stderr, "curiouserly, i got r = %d\n", r);
		FORI(5) fprintf(stderr, "\tt[%d] = %g\n", i, t[i]);
		return false;
	}
}

static bool load_camera_position_f(state *e, FILE *f)
{
	setlocale(LC_NUMERIC, "C"); // who ordered that?!
	float t[5] = {-42, -42, -42, -42, -42};
	int r = 0;
	FORI(5) r += fscanf(f, "%g", t + i);
	if (r == 5) {
		FORI(4) e->view_quat[i] = t[i];
		e->view_scale = t[4];
		FORI(3) e->view_quat_diff[i] = 0;
		e->view_quat_diff[3] = 1;
		expose_event(e->widget, NULL, e);
		return true;
	} else {
		fprintf(stderr, "strangely, i got r = %d\n", r);
		FORI(5) fprintf(stderr, "\tt[%d] = %g\n", i, t[i]);
		int f_a = feof(f);
		int f_b = ferror(f);
		fprintf(stderr, "feof = %d, ferror = %d\n", f_a, f_b);
		return false;
	}
}
static bool load_camera_position(state *e, char *fname)
{
	FILE *f = xfopen(fname, "r");
	bool r = load_camera_position_f(e, f);
	xfclose(f);
	return r;
}

static void save_lights_status_f(state *e, FILE *f)
{
	fprintf(f, "%d\n", MAX_LIGHTS);
	FORI(MAX_LIGHTS) {
		fprintf(f, "%d", e->light_active[i]);
		fprintf(f, "\n\t");
		FORJ(4) fprintf(f, " %g", e->light_ambient[i][j]);
		fprintf(f, "\n\t");
		FORJ(4) fprintf(f, " %g", e->light_diffuse[i][j]);
		fprintf(f, "\n\t");
		FORJ(4) fprintf(f, " %g", e->light_specular[i][j]);
		fprintf(f, "\n\t");
		FORJ(4) fprintf(f, " %g", e->light_position[i][j]);
		fprintf(f, "\n");
	}
}
static void save_lights_status(state *e, char *fname)
{ FILE *f = xfopen(fname, "w"); save_lights_status_f(e, f); xfclose(f); }

static bool load_lights_status_f(state *e, FILE *f)
{
	setlocale(LC_NUMERIC, "C"); // who ordered that?!
	int ml = 0, r; r=fscanf(f, "%d\n", &ml);
	if (ml > MAX_LIGHTS) return false;
	r = 0;
	FORI(MAX_LIGHTS) {
		int tmp = false;
		r += fscanf(f, "%d", &tmp);
		e->light_active[i] = tmp;
		FORJ(4) r += fscanf(f, "\n%g", e->light_ambient[i]+j);
		FORJ(4) r += fscanf(f, "\n%g", e->light_diffuse[i]+j);
		FORJ(4) r += fscanf(f, "\n%g", e->light_specular[i]+j);
		FORJ(4) r += fscanf(f, "\n%g", e->light_position[i]+j);
	}
	fprintf(stderr, "read %d nums (%d)\n", r, ml * (17));
	return true;
}
static bool load_lights_status(state *e, char *fname)
{
	FILE *f = xfopen(fname, "r");
	bool r = load_lights_status_f(e, f);
	xfclose(f);
	return r;
}

static void save_material_status_f(state *e, FILE *f)
{
	FORI(4) fprintf(f, "%g%c", e->material_ambient[i], i==3?'\n':' ');
	FORI(4) fprintf(f, "%g%c", e->material_diffuse[i], i==3?'\n':' ');
	FORI(4) fprintf(f, "%g%c", e->material_specular[i], i==3?'\n':' ');
	FORI(4) fprintf(f, "%g%c", e->material_emission[i], i==3?'\n':' ');
	fprintf(f, "%g\n", e->material_shininess);
}
static void save_material_status(state *e, char *fname)
{ FILE *f = xfopen(fname, "w"); save_material_status_f(e, f); xfclose(f); }

static bool load_material_status_f(state *e, FILE *f)
{
	setlocale(LC_NUMERIC, "C"); // who ordered that?!
	int r = 0;
	FORJ(4) r += fscanf(f, "%g\n", e->material_ambient + j);
	FORJ(4) r += fscanf(f, "%g\n", e->material_diffuse + j);
	FORJ(4) r += fscanf(f, "%g\n", e->material_specular + j);
	FORJ(4) r += fscanf(f, "%g\n", e->material_emission + j);
	r += fscanf(f, "%g\n", &e->material_shininess);
	fprintf(stderr, "read %d nums (%d)\n", r, 17);
	return true;
}
static bool load_material_status(state *e, char *fname)
{
	FILE *f = xfopen(fname, "r");
	bool r = load_material_status_f(e, f);
	xfclose(f);
	return r;
}

static void save_all_status(state *e, char *fname)
{
	FILE *f = xfopen(fname, "w");
	save_camera_position_f(e, f);
	save_material_status_f(e, f);
	save_lights_status_f(e, f);
	xfclose(f);
}

static bool load_all_status(state *e, char *fname)
{
	FILE *f = xfopen(fname, "r");
	bool r = true;
	r = r && load_camera_position_f(e, f);
	r = r && load_material_status_f(e, f);
	r = r && load_lights_status_f(e, f);
	xfclose(f);
	return r;
}

static int command_partbranch(state *e)
{
	// bureaucracy
	assert(e->cmdline_oldmode == MODE_MSTREE);

	tr_tree *t = e->current_ms_tree;
	tr_region *r = e->current_ms_region;
	fprintf(stderr, "fem la gràfica de la regió %td (val = %g) en amunt\n",
			r-e->current_ms_tree->the_regions, r->data.f);

	fprintf(stderr, "current_ms_tree = %d points, %d regions\n",
			t->number_of_things, t->number_of_regions);

	int cx = 0; while (r != t->root) { cx += 1; r = r->mother; }
	fprintf(stderr, "%d superfícies cap a l'arrel\n", cx);
	r = e->current_ms_region;
	int dx = 0; while (r) { dx += 1; r = r->daughter; }
	fprintf(stderr, "... i %d superfícies cap avall\n", dx);
	DXMALLOC(int, tidx, cx + dx + 1);

	r = e->current_ms_region;
	tidx[cx] = r - t->the_regions;
	r = e->current_ms_region;
	int k = cx;
	while (r != t->root) { r = r->mother; tidx[--k] = r-t->the_regions; }
	r = e->current_ms_region;
	assert(!k);
	k = cx;
	while(r) { r = r->daughter; tidx[++k] = r-t->the_regions; }
	assert(k==cx+dx);

	int nqpoints = cx + dx;
	float val[nqpoints];
	float che[nqpoints];
	int chemax = 0;
	FORI(nqpoints)
	{
		tr_region *s = t->the_regions + tidx[i];
		fprintf(stderr, "%dth%s: %d \"%g\" %d\n", i, i==cx?"(mine)":"",
				tidx[i],
				t->the_regions[tidx[i]].data.f,
				t->the_regions[tidx[i]].volume);
		val[i] = s->data.f;
		//che[i] = s->volume;
		che[i] = significativitat_desolneux(e->enfa, s->volume, s->data.f);
		if (che[i] > che[chemax]) chemax = i;
	}


	void raw_multiplot_bounds_s(int nplots,
			float xmin, float xmax,
			float ymin, float ymax,
			float **x, float **y, int *n, char *s[], char *t);
	float *xcoords[] = {val};
	float *ycoords[] = {che};
	int enn[] = {nqpoints};
	float max_che, min_che;
	float max_val, min_val;
	get_min_maxf(che, nqpoints, &min_che, &max_che);
	get_min_maxf(val, nqpoints, &min_val, &max_val);
	char *titles[] = {"significativitat desolneux..."};
	raw_multiplot_bounds_s(1, min_val, max_val, min_che, max_che,
			xcoords, ycoords, enn, titles, NULL);


	xfree(tidx);

	return tidx[chemax];
}

static int command_chebranch(state *e, int ridx, int nqpoints)
{
	// bureaucracy
	assert(ridx >= 0);
	assert(ridx < e->t->n);
	int bidx = -1;
	FORI(e->number_of_branches)
	{
		region *r = e->t->the_regions + e->branches_base[i];
		FORJ(e->branches_size[i])
		{
			if (r - e->t->the_regions ==  ridx)
			{ bidx = i; goto found; }
			r = r->parent;
		}
	}
found:	assert(bidx >= 0);
	assert(bidx == e->branches_which[ridx]);
	if (bidx != e->branches_which[ridx]) error("caca de la ben grossa");

	// fill branch indices
	int stb = e->branches_size[bidx];
	int thisbranch[stb];
	thisbranch[0] = e->branches_base[bidx];
	{
		region *r = e->t->the_regions + e->branches_base[bidx];
		FORI(stb - 1) {
			r = r->parent;
			thisbranch[i+1] = r - e->t->the_regions;
		}
	}

	// quantize this branch
	if (nqpoints < 6) nqpoints = 6;
	if (nqpoints > stb) nqpoints = stb;
	int qpoints[nqpoints];
	FORI(nqpoints) {
		int fi = (stb - 1) * i;
		fi /= nqpoints - 1;
		assert(fi >= 0);
		assert(fi < stb);
		qpoints[i] = thisbranch[fi];
		if (qpoints[i] < 0 || qpoints[i] >= e->t->n)
			error("bad qpoints[%d] == %d\n", i, qpoints[i]);
	}

	// data storage space
	float euclidean_volume[nqpoints];
	float euclidean_perimeter[nqpoints];
	float ponderated_perimeter[nqpoints];
	float val[nqpoints];
	float che[nqpoints];
	FORI(nqpoints) euclidean_volume[i] = NAN;
	FORI(nqpoints) euclidean_perimeter[i] = NAN;
	FORI(nqpoints) ponderated_perimeter[i] = NAN;
	FORI(nqpoints) val[i] = e->t->the_regions[qpoints[i]].value;
	FORI(nqpoints) che[i] = NAN;

	// process selected surfaces
	int minche = 0;
	FORI(nqpoints) {
		region *r = e->t->the_regions + qpoints[i];
		assert(r >= e->t->the_regions);
		assert(r < e->t->the_regions + e->t->n);
		if (e->t->root == r) continue;

		surface_triangulation st[1];
		find_cleant_triangulation_of_border_alpha(st, e->t, r,
				e->mmc_option, 0.5);
		fes_connectivity(st);
		orient_consistently(st);
		halfpixel_inline(st);
		float pinfty[] = {-1, -1, -1};
		euclidean_volume[i] = enclosed_volume(st, pinfty);
		euclidean_perimeter[i] = isotropic_riemannian_area(st, NULL);
		ponderated_perimeter[i] = isotropic_riemannian_area(st, e->g);
		free_surface_triangulation(st);
		che[i] = fabs(ponderated_perimeter[i]/euclidean_volume[i]);
		if (che[i] < che[minche]) minche = i;
	}
	minche = qpoints[minche];

	// draw things and cleanup
	void raw_multiplot_bounds_s(int nplots,
			float xmin, float xmax,
			float ymin, float ymax,
			float **x, float **y, int *n, char *s[], char *t);
	float *xcoords[] = {val};
	float *ycoords[] = {che};
	int enn[] = {nqpoints};
	float max_che, min_che;
	float max_val, min_val;
	get_min_maxf(che, nqpoints, &min_che, &max_che);
	get_min_maxf(val, nqpoints, &min_val, &max_val);
	char *titles[] = {"cheeger ratio..."};
	raw_multiplot_bounds_s(1, min_val, max_val, min_che, max_che,
			xcoords, ycoords, enn, titles, NULL);

	return minche;
}

static void process_light_command(state *e, char *s)
{
	int lindex, r, a, b;
	char cmd;
	float x[4];
	r = sscanf(s, "%d %c %g %g %g %g\n", &lindex, &cmd, x, x+1, x+2, x+3);
	fprintf(stderr, "processing ligth command \"%s\", r = %d\n", s, r);
	if (r != 2 && r != 6 && r != 3) return;
	if (r == 2 && (cmd != '0' && cmd != '1')) return;
	if (r == 3 && (cmd != 'v')) return;
	if (lindex < -1 || lindex >= MAX_LIGHTS) return;
	if(lindex==-1){a=0;b=MAX_LIGHTS;}else{a=lindex;b=a+1;}
	for (int i = a; i < b; i++) {//ugly hack
		switch(cmd) {
		case '0': e->light_active[i] = false; break;
		case '1': e->light_active[i] = true; break;
		case 'a': FORJ(4) e->light_ambient[i][j] = x[j]; break;
		case 'd': FORJ(4) e->light_diffuse[i][j] = x[j]; break;
		case 's': FORJ(4) e->light_specular[i][j] = x[j]; break;
		case 'p': FORJ(4) e->light_position[i][j] = x[j]; break;
		case 'm': FORJ(4) e->light_model_ambient[j] = x[j]; break;
		case 'v': e->light_local_view[0] = *x; break;
		default: return;
		}
	}
}

static void process_material_command(state *e, char *s)
{
	char cmd;
	float x[4];
	int r = sscanf(s, "%c %g %g %g %g\n", &cmd, x, x+1, x+2, x+3);
	fprintf(stderr, "processing material command \"%s\", r = %d\n", s, r);
	if (r != 2 && r != 5) return;
	if (r == 2 && (cmd != 'h')) return;
	switch(cmd) {
	case 'a': FORJ(4) e->material_ambient[j] = x[j]; break;
	case 'd': FORJ(4) e->material_diffuse[j] = x[j]; break;
	case 's': FORJ(4) e->material_specular[j] = x[j]; break;
	case 'e': FORJ(4) e->material_emission[j] = x[j]; break;
	case 'h': e->material_shininess = *x; break;
	default: return;
	}
}

static void process_clip_command(state *e, char *s)
{
	char cmd;
	float x[4];
	int r = sscanf(s, "%g %g %g %g\n", x, x+1, x+2, x+3);
	fprintf(stderr, "processing clip command \"%s\", r = %d\n", s, r);
	//if (r != 1 && r != 4) return;
	if (r == 1)
		e->clipplane_active[0] = *x > 0;
	if (r == 4)
		FORI(4)
			e->clipplane_equation[0][i] = x[i];
	//switch(cmd) {
	//case 'a': FORJ(4) e->material_ambient[j] = x[j]; break;
	//case 'd': FORJ(4) e->material_diffuse[j] = x[j]; break;
	//case 's': FORJ(4) e->material_specular[j] = x[j]; break;
	//case 'e': FORJ(4) e->material_emission[j] = x[j]; break;
	//case 'h': e->material_shininess = *x; break;
	//default: return;
	//}
}


/******************************
 * END TABLE OF KEY FUNCTIONS *
 ******************************/

// global_static_table_of_key_functions {{{1

static void (* const global_static_table_of_key_functions[][0x100])(state *) =
{
#define G_S_T_KF_COMMON ['q'] = kf_quit,\
			['o'] = kf_reverse_orientation,\
			['O'] = kf_toggle_surface_outlines,\
			['v'] = kf_toggle_selection_visibility,\
			['T'] = kf_toggle_tree_visibility,\
			['N'] = kf_rotate_point_visibility,\
			['D'] = kf_toggle_dotty_style,\
			['s'] = kf_toggle_wireframe,\
			['w'] = kf_use_wireframe,\
			['y'] = kf_yank_active_surface,\
			['Z'] = kf_save_selection_coh,\
			['x'] = kf_delete_last_selection,\
			['g'] = kf_increase_gamma,\
			['G'] = kf_decrease_gamma,\
			['V'] = kf_compute_nice_gradient,\
			['H'] = kf_rotate_hackmode,\
			['B'] = kf_toggle_white_background,\
			['P'] = kf_toggle_stats,\
			['E'] = kf_toggle_planes,\
			['R'] = kf_rotate_projection_type,\
			['4'] = kf_randomize_materialcolors,\
			['t'] = kf_mode_tos,\
			['m'] = kf_mode_mstree,\
			['p'] = kf_mode_point,\
			['a'] = kf_mode_malist,\
			[':'] = kf_mode_cmdline
	[MODE_TOS] = { G_S_T_KF_COMMON,
		['r'] = kf_tos_root,
		['j'] = kf_tos_daughter,
		['k'] = kf_tos_mother,
		['l'] = kf_tos_sister,
		['h'] = kf_tos_sister_back,
		['u'] = kf_tos_daughter_fast,
		['J'] = kf_tos_bottom_of_branch,
		['K'] = kf_tos_top_of_branch,
		['X'] = kf_tos_delete_current_shape,
		['Y'] = kf_tos_delete_current_shape_downwards,
		['W'] = kf_tos_save_tc,\
		['i'] = kf_tos_mother_fast,
		['0'] = kf_tos_option_zero,
		['1'] = kf_tos_option_one,
		['2'] = kf_tos_option_two
	},
	[MODE_MSTREE] = { G_S_T_KF_COMMON,
		['r'] = kf_mstree_root,
		['j'] = kf_mstree_daughter,
		['k'] = kf_mstree_mother,
		['l'] = kf_mstree_sister,
		['h'] = kf_mstree_sister_back
	},
	[MODE_POINT] = { G_S_T_KF_COMMON,
		['k'] = kf_point_dir_ym,
		['j'] = kf_point_dir_yp,
		['h'] = kf_point_dir_xm,
		['l'] = kf_point_dir_xp,
		['d'] = kf_point_dir_zm,
		['u'] = kf_point_dir_zp,
		[' '] = kf_point_put,
		['\b'] = kf_point_unput,
		['K'] = kf_point_make_center,
		['W'] = kf_point_save_asc,
	},
	[MODE_MALIST] = { G_S_T_KF_COMMON,
	},
	[MODE_CMDLINE] = {
		// THIS MODE IS PROCESSED SEPARATELY BY "cmdline_next_key"
	}
#undef G_S_T_KF_COMMON
};

static void check_stof(void)
{
#define t global_static_table_of_key_functions
	printf("checking stof:\n");
	FORI((int)LENGTH(t))
		FORJ((int)LENGTH(t[0]))
			if (t[i][j])
				printf("\tentry [%d][%d] ('%c') \"%p\"\n",
						i, j, j, t[i][j]);
#undef t
}


// command line management {{{1

static void cmdline_update(state *e, char *fmt, ...)
{
	if (fmt)
	{
		char buf[CMDLINE_LENGTH]; buf[0] = ' ';
		va_list argp;
		va_start(argp, fmt);
		int nbuf = vsnprintf(buf+1, CMDLINE_LENGTH-2, fmt, argp);
		va_end(argp);
		fprintf(stderr, "STATUS LINE MESSAGE: \"%s\"\n", buf+1);
		memcpy(e->cmdline_string, buf, nbuf+1);
		e->cmdline_string[nbuf+1] = 0;
	}
#ifdef USE_STATUSBAR_WIDGET
	gtk_statusbar_push(GTK_STATUSBAR(e->cmdline_w),
			1, e->cmdline_string);
#else
	gtk_label_set_text(GTK_LABEL(e->cmdline_w),
			e->cmdline_string);
#endif
}

// embryonic support for a vi-like ":" command line
// so far, only ":q" and a few options are implemented
// TODO : write a parser instead of hacking everything in
static void cmdline_process(state *e)
{
	assert(e->cmdline_pos < CMDLINE_LENGTH);
	assert(!e->cmdline_string[e->cmdline_pos]);
	fprintf(stderr, "processing command line \"%s\"\n", e->cmdline_string);
	char *s = e->cmdline_string;
	if (0 == strcmp(s, ":qa")) { // quit badly
		exit(EXIT_FAILURE);
	} else if (0 == strcmp(s, ":q") || 0 == strcmp(s, ":exit")) { // quit
		kf_quit(e);
	} else if (0 == strncmp(s, ":sleep ", 7)) { // sleep
		int d = atoi(s+7);
		cmdline_update(e, "sleeping for %d seconds...", d);
		idle(e->widget);
		sleep(d);
		cmdline_update(e, "done!");
	} else if (0 == strncmp(s, ":stride ", 8)) { // set unreal_five
		int n = atoi(s+8);
		if (n > 0)
		{
			e->unreal_five = n;
			cmdline_update(e, "new stride = %d", n);
		}
	} else if (0 == strncmp(s, ":ridx ", 6)) { // select region manually
		int nridx = atoi(s+6);
		if (nridx < 0 || nridx >= e->t->n)
			cmdline_update(e, "ERROR: bad region: %d", s+6);
		else {
			e->current_region = e->t->the_regions + nridx;
			treat_new_tos_region(e);
		}
	} else if (0 == strncmp(s, ":branch", 7)) { // display branch
		int ridx = e->current_region - e->t->the_regions;
		int bidx = e->branches_which[ridx];
		int bb = e->branches_base[bidx];
		int bl = e->branches_size[bidx];
		cmdline_update(e, "region %d inside branch %d of size %d",
				ridx, bb, bl);
	} else if (0 == strncmp(s, ":chebranch ", 11)) { // graph cheegerity
		int ridx = e->current_region - e->t->the_regions;
		int nqpoints = atoi(s+11);
		cmdline_update(e, "checomp of %d points from region %d...",
							nqpoints, ridx);
		int minche = command_chebranch(e, ridx, nqpoints);
		assert(minche >= 0);
		assert(minche < e->t->n);
		e->current_region = e->t->the_regions + minche;
		treat_new_tos_region(e);
	} else if (0 == strcmp(s, ":partbranch")) { // branch significativity
		int chemax = command_partbranch(e);
		if (chemax >= 0) {
			assert(chemax < e->current_ms_tree->number_of_regions);
			e->current_ms_region = e->current_ms_tree->the_regions;
			e->current_ms_region += chemax;
			treat_new_mstree_region(e);
		}
	} else if (0 == strncmp(s, ":hide", 5)) { // hide status bar
		gtk_widget_hide(e->cmdline_w);
	} else if (0 == strncmp(s, ":resize ", 8)) {
		int w, h, r = sscanf(s+8, "%d %d", &w, &h);
		if (r == 2) gtk_window_resize(GTK_WINDOW(e->window), w, h);
	} else if (0 == strcmp(s, ":fullscreen")) {
		gtk_window_fullscreen(GTK_WINDOW(e->window));
		//configure_event(e->widget, NULL, e);
		expose_event(e->widget, NULL, e);
		gdk_flush();
		while (g_main_context_iteration(NULL, FALSE));
	} else if (0 == strcmp(s, ":unfullscreen")) {
		gtk_window_unfullscreen(GTK_WINDOW(e->window));
	} else if (0 == strcmp(s, ":maximize")) {
		gtk_window_maximize(GTK_WINDOW(e->window));
	} else if (0 == strcmp(s, ":unmaximize")) {
		gtk_window_unmaximize(GTK_WINDOW(e->window));
	} else if (0 == strcmp(s, ":undecorate")) {
		gtk_window_set_decorated(GTK_WINDOW(e->window), FALSE);
	} else if (0 == strcmp(s, ":redecorate")) {
		gtk_window_set_decorated(GTK_WINDOW(e->window), TRUE);
	} else if (0 == strncmp(s, ":h", 2)) { // display help
		cmdline_update(e, "don't panic!");
	} else if (0 == strncmp(s, ":B ", 3)) { // read "binary" image
		imatge B[1], *v = e->voxel_drawing_image;
		carrega_imatge_general(B, s+3);
		if (same_sizeP(B, v)) {
			memcpy(v->grayplane, B->grayplane,
					v->mida * sizeof(gris));
			kf_mode_point(e);
		} else {
			cmdline_update(e, "ERROR: bad image \"%s\"", s+3);
		}
		allibera_imatge(B);
	} else if (0 == strncmp(s, ":wp ", 4)) { // write "binary" image
		desa_imatge_general(e->voxel_drawing_image, s+4);
		cmdline_update(e, "wrote binary image to file \"%s\"", s+4);
	} else if (0 == strcmp(s, ":minshapeB")) { // smallest shape > B
		int mi = min_shape_containingB(e->t, e->voxel_drawing_image);
		assert(mi >= 0);
		assert(mi < e->t->n);
		e->current_region = e->t->the_regions + mi;
		treat_new_tos_region(e);
	} else if (0 == strncmp(s, ":savepos ", 9)) { // save 3D camera view
		save_camera_position(e, s+9);
		cmdline_update(e, "wrote camera position to file \"%s\"", s+9);
	} else if (0 == strncmp(s, ":loadpos ", 9)) { // load 3D camera view
		char *msg = load_camera_position(e, s+9) ?
			"read camera position from file \"%s\""
			: "ERROR: could not read camera pos from file \"%s\"";
		cmdline_update(e, msg , s+9);
	} else if (0 == strncmp(s, ":savelit ", 9)) { // save 3D camera view
		save_lights_status(e, s+9);
		cmdline_update(e, "wrote lights status to file \"%s\"", s+9);
	} else if (0 == strncmp(s, ":loadlit ", 9)) { // load 3D camera view
		char *msg = load_lights_status(e, s+9) ?
			"read ligths status from file \"%s\""
			: "ERROR: could not read lit status from file \"%s\"";
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
		cmdline_update(e, msg , s+9);
	} else if (0 == strncmp(s, ":savemat ", 9)) { // save 3D camera view
		save_material_status(e, s+9);
		cmdline_update(e, "wrote material status to file \"%s\"", s+9);
	} else if (0 == strncmp(s, ":loadmat ", 9)) { // load 3D camera view
		char *msg = load_material_status(e, s+9) ?
			"read material status from file \"%s\""
			: "ERROR: could not read mat status from file \"%s\"";
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
		cmdline_update(e, msg , s+9);
	} else if (0 == strncmp(s, ":saveall ", 9)) { // save 3D camera view
		save_all_status(e, s+9);
		cmdline_update(e, "wrote all status to file \"%s\"", s+9);
	} else if (0 == strncmp(s, ":loadall ", 9)) { // load 3D camera view
		char *msg = load_all_status(e, s+9) ?
			"read all status from file \"%s\""
			: "ERROR: could not read all status from file \"%s\"";
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
		cmdline_update(e, msg , s+9);
	} else if (0 == strncmp(s, ":loadsingle ", 12)) {
	} else if (0 == strncmp(s, ":sp", 3)) { // save 3D view (default file)
		save_camera_position(e, DEFAULT_CAMERA_POS_FNAME);
	} else if (0 == strncmp(s, ":lp", 3)) { // load 3D view (default file)
		load_camera_position(e, DEFAULT_CAMERA_POS_FNAME);
	} else if (0 == strncmp(s, ":setpos ", 8)) { // set 3D camera view
		set_camera_position(e, s+8);
	} else if (0 == strncmp(s, ":trackview", 6)) { // view trackball
		cmdline_update(e, "trackball [pos ; diff] = [%g %g ; %g %g]",
				e->trackball_begin_x, e->trackball_begin_y,
				e->trackball_dx, e->trackball_dy);
	} else if (0 == strncmp(s, ":quatview", 6)) { // view trackball
		cmdline_update(e, "[q ; qdiff] = [%g %g %g %g ; %g %g %g %g]",
				e->view_quat[0], e->view_quat[1],
				e->view_quat[2], e->view_quat[3],
				e->view_quat_diff[0], e->view_quat_diff[1],
				e->view_quat_diff[2], e->view_quat_diff[3]);
	} else if (0 == strncmp(s, ":l ", 3)) {
		process_light_command(e, s+3);
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":c ", 3)) {
		process_clip_command(e, s+3);
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":m ", 3)) {
		process_material_command(e, s+3);
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":expose", 7)) {
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":shot ", 6)) {
		save_screenshot(e, s + 6);
		cmdline_update(e, "wrote screenshot to file \"%s\"", s + 6);
	} else if (0 == strncmp(s, ":size", 5)) {
		int w = e->drawing_area->allocation.width;
		int h = e->drawing_area->allocation.height;
		cmdline_update(e, "sizes = %d %d", w, h);
	} else if (0 == strncmp(s, ":ilit", 5)) {
		FORI(MAX_LIGHTS) FORJ(3) e->light_position[i][j] *= -1;
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":drawnormals", 12)) {
		e->selection_normals_draw = !e->selection_normals_draw;
		if (e->selection_normals_draw && !e->selection->t->n)
			compute_selection_normals(e);
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":normals", 8)) {
		if (!e->selection->t->n)
			compute_selection_normals(e);
		else
			FORI(e->selection->number_of_pieces)
			{
				xfree(e->selection->t[i].n);
				e->selection->t[i].n = NULL;
		}
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":smotion ", 9)) {
		float lambda = atof(s+9);
		FORI(e->selection->number_of_pieces)
			surface_normal_motion(e->selection->t + i, lambda);
		e->selection_updated = true;
		dump_objects(e);
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":axes", 5)) {
		e->show_drawing_of_axes = ! e->show_drawing_of_axes;
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":ldots", 6)) {
		e->lights_displayed_as_dots = ! e->lights_displayed_as_dots;
		expose_event(e->widget, NULL, e);
	} else if (0 == strncmp(s, ":echo", 5)) { // echo
		cmdline_update(e, s + (' '!=s[5]?5:6));
	} else { // bad command
		cmdline_update(e, "ERROR: bad command: %s", s+1);
	}

	// TODO : write complex stuff here
}

static void cmdline_next_key(state *e, int key)
{
	assert(key >= 0);
	assert(key < 0x100);

	//fprintf(stderr, "cmdline char %d {%x} '%c'\n", key, key, key);

	switch(key) {
	case '\n':	cmdline_process(e);
	case '\e':	e->mode = e->cmdline_oldmode;
	case 0:		break;
	default:	if (e->cmdline_pos < CMDLINE_LENGTH - 2)
			{
				e->cmdline_string[e->cmdline_pos] = key;
				e->cmdline_pos += 2;
	case '\b':	    if (e->cmdline_pos > 0) e->cmdline_pos -= 1;
				e->cmdline_string[e->cmdline_pos] = 0;
				cmdline_update(e, NULL);
			} else assert(false);
	} // purposely spaghettified
}

// OpenGL debugging {{{1
static void
print_gl_config_attrib (GdkGLConfig *glconfig,
                        const gchar *attrib_str,
                        int          attrib,
                        gboolean     is_boolean)
{
	int value;

	g_print ("%s = ", attrib_str);
	if (gdk_gl_config_get_attrib (glconfig, attrib, &value))
	{
		if (is_boolean)
			g_print ("%s\n", value == TRUE ? "TRUE" : "FALSE");
		else
			g_print ("%d\n", value);
	}
	else
		g_print ("*** Cannot get %s attribute value\n", attrib_str);
}

static void
examine_gl_config_attrib (GdkGLConfig *glconfig)
{
	g_print ("\nOpenGL visual configurations :\n\n");

	g_print ("gdk_gl_config_is_rgba (glconfig) = %s\n",
		gdk_gl_config_is_rgba (glconfig) ? "TRUE" : "FALSE");
	g_print ("gdk_gl_config_is_double_buffered (glconfig) = %s\n",
		gdk_gl_config_is_double_buffered (glconfig) ? "TRUE" : "FALSE");
	g_print ("gdk_gl_config_is_stereo (glconfig) = %s\n",
		gdk_gl_config_is_stereo (glconfig) ? "TRUE" : "FALSE");
	g_print ("gdk_gl_config_has_alpha (glconfig) = %s\n",
		gdk_gl_config_has_alpha (glconfig) ? "TRUE" : "FALSE");
	g_print ("gdk_gl_config_has_depth_buffer (glconfig) = %s\n",
		gdk_gl_config_has_depth_buffer (glconfig) ? "TRUE" : "FALSE");
	g_print ("gdk_gl_config_has_stencil_buffer (glconfig) = %s\n",
		gdk_gl_config_has_stencil_buffer (glconfig) ? "TRUE" : "FALSE");
	g_print ("gdk_gl_config_has_accum_buffer (glconfig) = %s\n",
		gdk_gl_config_has_accum_buffer (glconfig) ? "TRUE" : "FALSE");

	g_print ("\n");
	print_gl_config_attrib (glconfig, SCL(GDK_GL_USE_GL),           TRUE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_BUFFER_SIZE),      FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_LEVEL),            FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_RGBA),             TRUE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_DOUBLEBUFFER),     TRUE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_STEREO),           TRUE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_AUX_BUFFERS),      FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_RED_SIZE),         FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_GREEN_SIZE),       FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_BLUE_SIZE),        FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_ALPHA_SIZE),       FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_DEPTH_SIZE),       FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_STENCIL_SIZE),     FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_ACCUM_RED_SIZE),   FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_ACCUM_GREEN_SIZE), FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_ACCUM_BLUE_SIZE),  FALSE);
	print_gl_config_attrib (glconfig, SCL(GDK_GL_ACCUM_ALPHA_SIZE), FALSE);

	g_print ("\n");
}

// GTK and OpenGL initialization {{{1
static void setup_gtk_and_gl_environments(state *e, int argc, char *argv[])
{
	GdkGLConfig *glconfig;
	gint major, minor;

	GtkWidget *window;
	GtkWidget *vbox;
	GtkWidget *drawing_area;
	GtkWidget *menu;
	GtkWidget *button;

	/* Initialize GTK. */
	gtk_init (&argc, &argv);

	/* Initialize GtkGLExt. */
	gtk_gl_init (&argc, &argv);


	/*
	 * Query OpenGL extension version.
	 */

	gdk_gl_query_version (&major, &minor);
	g_print ("\nOpenGL extension version - %d.%d\n",
			major, minor);

	/*
	 * Configure OpenGL-capable visual.
	 */

	/* Try double-buffered visual */
	glconfig = gdk_gl_config_new_by_mode ( (GdkGLConfigMode)
			(GDK_GL_MODE_RGB    |
			 GDK_GL_MODE_DEPTH  |
			 GDK_GL_MODE_ALPHA  |
			 GDK_GL_MODE_DOUBLE));
	glconfig = NULL;
	if (glconfig == NULL)
	{
		g_print ("*** Cannot find the double-buffered visual.\n");
		g_print ("*** Trying single-buffered visual.\n");

		/* Try single-buffered visual */
		glconfig = gdk_gl_config_new_by_mode ( (GdkGLConfigMode)
				(GDK_GL_MODE_RGB   |
				 GDK_GL_MODE_DEPTH));
		if (glconfig == NULL)
		{
			g_print ("*** No OpenGL-capable visual found.\n");
			exit (1);
		}
	}

	examine_gl_config_attrib (glconfig);

	/*
	 * Top-level window.
	 */

	window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	e->window = window;
	gtk_window_set_title (GTK_WINDOW (window), "glshapest");

	/* Get automatically redrawn if any of their children moves. */
	gtk_container_set_reallocate_redraws (GTK_CONTAINER (window), TRUE);

	g_signal_connect (G_OBJECT (window), "delete_event",
			G_CALLBACK (gtk_main_quit), NULL);

	/*
	 * VBox.
	 */

	vbox = gtk_vbox_new (FALSE, 0);
	gtk_container_add (GTK_CONTAINER (window), vbox);
	gtk_widget_show (vbox);

	/*
	 * Drawing area for drawing OpenGL scene.
	 */

	drawing_area = gtk_drawing_area_new ();
	e->drawing_area = drawing_area;
	gtk_widget_set_size_request (drawing_area, 600, 600);

	/* Set OpenGL-capability to the widget. */
	gtk_widget_set_gl_capability (drawing_area,
			glconfig,
			NULL,
			TRUE,
			GDK_GL_RGBA_TYPE);

	gtk_widget_add_events (drawing_area,
			GDK_BUTTON1_MOTION_MASK    |
			GDK_BUTTON2_MOTION_MASK    |
			GDK_BUTTON_PRESS_MASK      |
			GDK_BUTTON_RELEASE_MASK    |
			GDK_VISIBILITY_NOTIFY_MASK);

	g_signal_connect_after (G_OBJECT (drawing_area), "realize",
			G_CALLBACK (realize), e);
	g_signal_connect (G_OBJECT (drawing_area), "configure_event",
			G_CALLBACK (configure_event), NULL);
	g_signal_connect (G_OBJECT (drawing_area), "expose_event",
			G_CALLBACK (expose_event), e);

	g_signal_connect (G_OBJECT (drawing_area), "button_press_event",
			G_CALLBACK (button_press_event), e);
	g_signal_connect (G_OBJECT (drawing_area), "button_release_event",
			G_CALLBACK (button_release_event), e);
	g_signal_connect (G_OBJECT (drawing_area), "motion_notify_event",
			G_CALLBACK (motion_notify_event), e);

	g_signal_connect (G_OBJECT (drawing_area), "map_event",
			G_CALLBACK (map_event), e);
	g_signal_connect (G_OBJECT (drawing_area), "unmap_event",
			G_CALLBACK (unmap_event), NULL);
	g_signal_connect (G_OBJECT (drawing_area), "visibility_notify_event",
			G_CALLBACK (visibility_notify_event), e);

	// XXX: swapping here kills the program!!!
	//g_signal_connect_swapped (G_OBJECT (window), "key_press_event",
	//		G_CALLBACK (key_press_event), e/*drawing_area*/);

	g_signal_connect(G_OBJECT (window), "key_press_event",
			G_CALLBACK (key_press_event), e);

	gtk_box_pack_start (GTK_BOX (vbox), drawing_area, TRUE, TRUE, 0);

	gtk_widget_show (drawing_area);

	/*
	 * coco's status line
	 */
#ifdef USE_STATUSBAR_WIDGET
	e->cmdline_w = gtk_statusbar_new();
	gtk_statusbar_push(GTK_STATUSBAR(e->cmdline_w),1,"bon dia!");
#else
	e->cmdline_w = gtk_label_new(" bon dia!");
	g_object_set(e->cmdline_w,
			"single-line-mode", TRUE,
			"xalign", 0.0,
			NULL);
#endif
	gtk_box_pack_start(GTK_BOX(vbox), e->cmdline_w, FALSE, FALSE, 2);
	//gtk_widget_show(e->cmdline_w);
	//gtk_widget_hide(e->cmdline_w);


	/*
	 * Simple quit button.
	 */

	//button = gtk_button_new_with_label ("Quit");

	//g_signal_connect (G_OBJECT (button), "clicked",
	//		G_CALLBACK (gtk_main_quit), NULL);

	//gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

	//gtk_widget_show (button);

#ifdef USE_FANCY_CMDLINE
	GtkWidget *fancy_cmdline = vte_terminal_new();
	vte_terminal_set_size(VTE_TERMINAL(fancy_cmdline), 80, 3);
	vte_terminal_set_color_background(VTE_TERMINAL(fancy_cmdline),
			&(GdkColor){255, 255, 255, 255});
	vte_terminal_feed(VTE_TERMINAL(fancy_cmdline), "bon dia!", 8);
	gtk_box_pack_start(GTK_BOX(vbox), fancy_cmdline, FALSE, FALSE, 2);
	gtk_widget_show(fancy_cmdline);
#endif //USE_FANCY_CMDLINE


	/*
	 * Show window.
	 */

	gtk_widget_show (window);
}

// argument handling {{{1

static region *find_starting_region(tree_of_regions *t)
{
	assert(t);
	assert(t->root);
	assert(t->root->first_child);
	region *r = t->root->first_child;
	if (t->mida > 0x1000)
	{
		int reasonable_volume = t->mida / 20;
		while (r->volume > reasonable_volume && r->first_child)
			r = r->first_child;
	}
	return r;
}

// the following function is to be seriously de-spaghettified:
static void load_arguments_into_state(state *e, int argc, char *argv[])
{
	// FIXME (does not seem to work: only accepts defaults)
	assert(argc == 5);
	tree_of_regions *t = e->tstatic;
	imatge *g = e->gstatic;
	printf("reading tc from \"%s\"...\n", argv[1]);
	carrega_tor_complet(t, argv[1]);
	printf("...read!\n");
	emplena_voxels(t);
	troba_un_pinf(t);
	float s = BAD_MAX(t->w, BAD_MAX(t->h, t->d)) - 1;
	int startingridx = atoi(argv[4]);
	printf("STATE pointer = %p\n", (void *)e);
	e->global_scale = 2.0/s;
	e->global_displacement[0] = - (t->w - 1) * e->global_scale / 2.0;
	e->global_displacement[1] = - (t->h - 1) * e->global_scale / 2.0;
	e->global_displacement[2] = - (t->d - 1) * e->global_scale / 2.0;
	reorder_children_by_volume_recursively(t, t->root);
	g->mida = 0;
	inicialitza_himatge(global_metric,t);FORI(t->mida)global_metric->grayplane[i]=1;
	if (0 == strcmp(argv[2], "."))
	{
		fprintf(stderr, "colors a lo cutre\n");
		imatge x;
		reconstruct(&x, t);
		fes_imatge_de_gradients_centrats(g, &x,
				perllonga_constant_lleig);
		allibera_imatge(&x);
	} else {
		fprintf(stderr, "colors presos de \"%s\"\n", argv[2]);
		carrega_imatge_general(g, argv[2]);
		FORI(t->mida)global_metric->grayplane[i]=g->grayplane[i];
	}
	if (g->mida != t->mida)
		error("size mismatch between tree and image");
	float min, max;
	get_min_maxf(g->grayplane, g->mida, &min, &max);
	e->imgrang = max;
	if (0 == strcmp(argv[3], "."))
	{
		init_state(e, t, g, NULL);
	} else {
		load_collection_of_hypersurfaces(e->hstatic, argv[3]);
		init_state(e, t, g, e->hstatic);
		e->tree_visibility = false;
	}
	//else init_state(e, t, NULL, NULL);


	if (startingridx >= 0)
		e->current_region = e->t->the_regions + startingridx;
	else
		e->current_region = find_starting_region(e->t);
	e->current_boundary = pick_boundary_triangulation(e, e->current_region);
	e->active_surface = e->current_boundary;
	//produce_triangulation();
}

static char *hyphen_c_hack = NULL;

// main program {{{1
static void parse_options(char *t[], int c, char *v[])
{
	int opt;
	extern char *optarg;
	while ((opt = getopt(c, v, "b:s:r:c:")) != -1) {
		switch (opt) {
		//case 'c':
		//case 'T':
		case 'b': t[2] = optarg; break;
		case 's': t[3] = optarg; break;
		case 'r': t[4] = optarg; break;
		case 'c': hyphen_c_hack = optarg; break;
		default: error("janderror (opt=%d)", opt);
		}
	}
}

static int main_full(int c, char *v[])
{
	char *t[5] = {"", "-", ".", ".", "-1"};
	parse_options(t, c, v);
	FORI(5) fprintf(stderr, "arg[%d]: \"%s\"\n", i, t[i]);
	state e[1];
	e->the_hyphen_c_argument = hyphen_c_hack; // hack


	setlocale(LC_NUMERIC, "C");

	load_arguments_into_state(e, 5, t);
	setup_gtk_and_gl_environments(e, c, v);

	/*
	 * Main loop.
	 */

	setlocale(LC_NUMERIC, "C"); // yes! it is necessary again
	gtk_main ();

	printf("sizeof(state) = %zu\n", sizeof e);
	printf("I'm outside of gtk_main!\n");
	clean_state(e);

	//check_stof();
	return 0;
}

int main(int c, char *v[])
{
	return ((c > 1) && helpargumentP(v[1])) ?
		fprintf(stderr, "usage:\n\t%s "
			"<in.tc [-b bg.asc] [-s pi.coh] [-r sidx] [-c cmd]\n",
			*v)
		: main_full(c, v);
}

// compilation line:
// gcc -std=gnu99 -Wall -Wextra -Wno-unused -g -rdynamic `pkg-config gtkglext-x11-1.0 --cflags --libs` old_glshapest.c gui/glshapest/trackball.o -I gui/glshapest -L.  -lqnm -lfftw3 -lgsl -lgslcblas -lm -o old_glshapest

// }}}1

// vim:set foldmethod=marker:
