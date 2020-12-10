#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "polyline.h"
#include "trackball.h"
#include "tmatrix.h"

Polyhedron* poly;

/*scene related variables*/
const float zoomspeed = 0.9;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 1.0;
int win_width = 800;
int win_height = 800;
float aspectRatio = win_width / win_height;
/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

/*IBFV related variables*/
//https://www.win.tue.nl/~vanwijk/ibfv/
#define	NPN 64
#define SCALE 4.0
int    Npat = 32;
int    iframe = 0;
float  tmax = win_width / (SCALE*NPN);
float  dmax = SCALE / win_width;
unsigned char *pixels;

#define DM  ((float) (1.0/(100-1.0)))

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void makePatterns(void);

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);




//*************** Final Project *************** 
#include "Streamline.h"
#include "StreamlinePlacementAlg.cpp"


float min_x = FLT_MAX, min_y = FLT_MAX, max_x = FLT_MIN, max_y = FLT_MIN;
bool streamlines = false, singularities = false;
bool foundStreamlines = false, foundSingleStreamline = false,  foundSingularities = false;

PolyLine Streamlines;
PolyLine SingleStreamline;
std::vector<Vertex*> Singularities;

void drawLineRecursive(Vertex* v, bool forward, int, bool);
int recursiveDepth = 200;
void findSingularities();
void extractSingularitiesQuad(Quad*);
void classifySingularity(Quad*, Vertex*);
void drawSingularities();
Vertex* getVertexAt(float xpos, float ypos);
Vertex* RKGetNextVertex(Vertex* initial, bool forward);



//*************************************** 





/*display vis results*/
void display_polyhedron(Polyhedron* poly);

/*display utilities*/

/*
draw a sphere
x, y, z are the coordiate of the dot
radius of the sphere 
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawDot(double x, double y, double z, double radius = 0.1, double R = 1.0, double G = 0.0, double B = 0.0) {

	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	GLfloat mat_diffuse[4];

	{
		mat_diffuse[0] = R;
		mat_diffuse[1] = G;
		mat_diffuse[2] = B;
		mat_diffuse[3] = 1.0;
	}

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	GLUquadric* quad = gluNewQuadric();

	glPushMatrix();
	glTranslatef(x, y, z);
	gluSphere(quad, radius, 50, 50);
	glPopMatrix();

	gluDeleteQuadric(quad);
}

/*
draw a line segment
width: the width of the line, should bigger than 0
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawLineSegment(LineSegment ls, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {

	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);
	glVertex3f(ls.start.x, ls.start.y, ls.start.z);
	glVertex3f(ls.end.x, ls.end.y, ls.end.z);
	glEnd();

	glDisable(GL_BLEND);
}

/*
draw a polyline
width: the width of the line, should bigger than 0
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawPolyline(PolyLine pl, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {
	
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);

	for (int i = 0; i < pl.size(); i++) {
		glVertex3f(pl[i].start.x, pl[i].start.y, pl[i].start.z);
		glVertex3f(pl[i].end.x, pl[i].end.y, pl[i].end.z);
	}

	glEnd();

	glDisable(GL_BLEND);
}

/******************************************************************************
Main program.
******************************************************************************/
int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../quadmesh_2D/new_vector_data/v1.ply", "r");
	//FILE* this_file = fopen("../quadmesh_2D/scalar_data/sin_function.ply", "r");

	poly = new Polyhedron(this_file);
	fclose(this_file);
	
	//find the max and min values of all scalars in the data set
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* v = poly->vlist[i];
		if (v->x < min_x)
			min_x = v->x;
		if (v->x > max_x)
			max_x = v->x;
		if (v->y < min_y)
			min_y = v->y;
		if (v->y > max_y)
			max_y = v->y;
	}


	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();

	/*prepare the noise texture for IBFV*/
	makePatterns();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	getVertexAt(4.9, 4.8);

	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}


/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor* zoom, radius_factor* zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor* zoom / aspectRatio, radius_factor* zoom / aspectRatio, 0.1, 1000);
	}


	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}


/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}


/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	//glDisable(GL_LIGHTING);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		{
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			//glColor3f(0, 0, 0);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Vertex* temp_v = this_poly->vlist[i];

		{
			GLUquadric* quad = gluNewQuadric();

			glPushMatrix();
			glTranslatef(temp_v->x, temp_v->y, temp_v->z);
			glColor4f(0, 0, 1, 1.0);
			gluSphere(quad, this_poly->radius * 0.01, 50, 50);
			glPopMatrix();

			gluDeleteQuadric(quad);
		}
	}
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	/*
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glDisable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 1.0);
		glVertex3d(temp_v->x, temp_v->y, 0.0);
	}
	glEnd();
	*/
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];
	Vertex* temp_v = temp_q->verts[0];

	drawDot(temp_v->x, temp_v->y, temp_v->z, this_poly->radius * 0.01, 1.0, 0.0, 0.0);

	//draw the streamline using the recursive call
	// draw both forward and backward
	
	if (foundSingleStreamline == false) {
		drawLineRecursive(temp_v, true, recursiveDepth, true);
		drawLineRecursive(temp_v, false, recursiveDepth, true);
		foundSingleStreamline = true;		SingleStreamline.clear();
		drawLineRecursive(temp_v, true, recursiveDepth, true);
		drawLineRecursive(temp_v, false, recursiveDepth, true);
		foundSingleStreamline = true;
	}
	


}

void drawAllStreamlines() {
	int subdivisions = 5;
	for (int i = 0; i < subdivisions + 1; i++) {
		for (int j = 0; j < subdivisions + 1; j++) {
			Vertex* v = getVertexAt(min_x + (float)i / (float)subdivisions * (max_x - min_x), min_y + (float)j / (float)subdivisions * (max_y - min_y));
			drawLineRecursive(v, true, recursiveDepth, false);
			drawLineRecursive(v, false, recursiveDepth, false);
		}
	}
	foundStreamlines = true;

}


void drawLineRecursive(Vertex* v, bool forward, int depthRemaining, bool single) {

	if (depthRemaining <= 0)
		return;
	Vertex* next_v;

	//Get the estimate of the next vertex from the Runge-Kutta algorithm
	if (forward) 
		next_v = RKGetNextVertex(v, true);
	else 
		next_v = RKGetNextVertex(v, false);
	
	//if the next point would go out of the data set, return 
	if (next_v == nullptr)
		return;

	//add the line segment to the streamline set
	LineSegment newLine(v->x, v->y, 0, next_v->x, next_v->y, 0);
	if (single)
		SingleStreamline.push_back(newLine);
	else
		Streamlines.push_back(newLine);

	//continue the line from the next point
	drawLineRecursive(next_v, forward, --depthRemaining, single);

}


Vertex* RKGetNextVertex(Vertex* initial, bool forward) {

	if ((float)initial->vx == 0 && (float)initial->vy == 0)
		return nullptr;

	//do Runge Kutta Method
	//time step will be found by using the magnitude of the vector at the point
	//ignore t? just use the x and y, increase them both as you go through the R-K method
	Vertex* v;
	
	int direction;
	direction = (forward ? 1 : -1);								//if forward = true, use the normal vectors, otherwise reverse them
	double timestep = .2 * direction , timeStepProportion;		//I think this isn't enough to make the forward/back work right, revisit later
	double nextx, nexty;

	
	icVector2 vec1(initial->vx, initial->vy);
	normalize(vec1);

	nextx = initial->x + timestep / 2.0 * vec1.x;
	nexty = initial->y + timestep / 2.0 * vec1.y;
	v = getVertexAt(nextx, nexty);
	if (v == nullptr)
		return nullptr;
	icVector2 vec2(v->vx, v->vy);
	normalize(vec2);

	nextx = initial->x + timestep / 2.0 * vec2.x;
	nexty = initial->y + timestep / 2.0 * vec2.y;
	v = getVertexAt(nextx, nexty);
	if (v == nullptr)
		return nullptr;
	icVector2 vec3(v->vx, v->vy);
	normalize(vec3);

	nextx = initial->x + timestep * vec3.x;
	nexty = initial->y + timestep * vec3.y;
	v = getVertexAt(nextx, nexty);
	if (v == nullptr)
		return nullptr;
	icVector2 vec4(v->vx, v->vy);
	normalize(vec4);

	Vertex* next_v;
	float finaly = initial->y + 1.0 / 6.0 * timestep * (vec1.y + 2.0 * vec2.y + 2.0 * vec3.y + vec4.y);
	float finalx = initial->x + 1.0 / 6.0 * timestep * (vec1.x + 2.0 * vec2.x + 2.0 * vec3.x + vec4.x);
	//printf("This x: %f This y: %f \nNext x: %f Next y: %f \n\n", initial->x, initial->y, finalx, finaly);
	if ((float)initial->x == finalx && (float)initial->y == finaly)
		return nullptr;
	return getVertexAt(finalx, finaly);

}

void findSingularities() {

	for (int i = 0; i < poly->nquads; i++) {
		Quad* this_q = poly->qlist[i];
		extractSingularitiesQuad(this_q);
	}

}

void extractSingularitiesQuad(Quad* q) {

	/*
	float max_x, max_y, min_x, min_y;
	max_x = max_y = FLT_MIN;
	min_x = min_y = FLT_MAX;
	*/

	double x1, x2, y1, y2;
	x1 = q->verts[2]->x;
	y1 = q->verts[2]->y;
	x2 = q->verts[0]->x;
	y2 = q->verts[0]->y;

	
	int k = 2;
	double fx1y1 = q->verts[k]->vx; //f(x, y) is the vx
	double fx2y1 = q->verts[(k + 1) % 4]->vx;
	double fx2y2 = q->verts[(k + 2) % 4]->vx;
	double fx1y2 = q->verts[(k + 3) % 4]->vx;
	double gx1y1 = q->verts[k]->vy; //g(x, y) is the vy
	double gx2y1 = q->verts[(k + 1) % 4]->vy;
	double gx2y2 = q->verts[(k + 2) % 4]->vy;
	double gx1y2 = q->verts[(k + 3) % 4]->vy;

	double a00 = fx1y1; 
	double a10 = fx2y1 - fx1y1; 
	double a01 = fx1y2 - fx1y1;
	double a11 = fx1y1 - fx2y1 - fx1y2 + fx2y2; 
	
	double b00 = gx1y1; 
	double b10 = gx2y1 - gx1y1; 
	double b01 = gx1y2 - gx1y1; 
	double b11 = gx1y1 - gx2y1 - gx1y2 + gx2y2; 
	
	double c00 = a11 * b00 - a00 * b11; 
	double c10 = a11 * b10 - a10 * b11; 
	double c01 = a11 * b01 - a01 * b11;


	double eqA = (-a11 * c10);
	double eqB = (-a11 * c00 - a01 * c10 + a10 * c01);
	double eqC = (a00 * c01 - a01 * c00);
	double discriminant = pow(eqB, 2) - (4 * eqA * eqC);
	//if (discriminant > 0)

	double s1 = (-eqB + sqrt(discriminant)) / (2 * eqA);
	double t1 = -(c00 / c01) - (c10 / c01) * s1;

	double s2 = (-eqB - sqrt(discriminant)) / (2 * eqA);
	double t2 = -(c00 / c01) - (c10 / c01 ) * s2;

	double finalx, finaly;
	bool valid = false;
	if(s2 > 0 && s2 < 1 || t2 > 0 && t2 < 1 &&
		s1 > 0 && s1 < 1 || t1 > 0 && t1 < 1)
		//printf("\n\ns1: %f, t1: %f, \n s2: %f, t2: %f", s1, t1, s2, t2);
	
	if (s1 > 0 && s1 < 1 && t1 > 0 && t1 < 1) {
		finalx = x1 + s1;
		finaly = y1 + t1;
		valid = true;
	}
	else if (s2 > 0 && s2 < 1 && t2 > 0 && t2 < 1) {
		finalx = x1 + s2;
		finaly = y1 + t2;
		valid = true;
	}
	
	if (valid) {
		Vertex * sing = new Vertex(finalx, finaly, 0);
		Singularities.push_back(sing);
		classifySingularity(q, sing);
	}
}

void classifySingularity(Quad* q, Vertex* v) {

	//walk around the vectors of the quad
	//as we walk, record the angle change between the vectors
	//if total angle change == 720
	//	singularity is a saddle
	//if 360
	//	if vectors pointing in
	//		singularity is a sink
	//	if vectors pointing out
	//		singularity is a source



	float totalAngle = 0;
	double theta, thisAngleChange, thisDir, nextDir;
	icVector2 thisVector, nextVector;
	for (int i = 0; i < 4; i++) {
		thisVector = icVector2(q->verts[i]->vx, q->verts[i]->vy);
		thisDir = atan(thisVector.y / thisVector.x);
		if (thisVector.x < 0)
			thisDir += PI;
		nextVector = icVector2(q->verts[(i + 1) % 4]->vx, q->verts[(i + 1) % 4]->vy);
		nextDir = atan(nextVector.y / nextVector.x);
		if (nextVector.x < 0)
			nextDir += PI;
		//thisMag = sqrt(pow(thisVector.x, 2) + pow(thisVector.y, 2));
		//nextMag = sqrt(pow(nextVector.x, 2) + pow(nextVector.y, 2));
		
		//theta = dot(thisVector, nextVector) / (thisMag * nextMag);
		//thisAngle = acos(theta);

		thisAngleChange = nextDir - thisDir;
		if (thisAngleChange > PI)
			thisAngleChange = -(2* PI - thisAngleChange);
		//if (thisAngleChange < -PI)
			//thisAngleChange = nextDir + thisDir;

		totalAngle += thisAngleChange;
		//printf("\nx: %f, y: %f - angle = %f", thisVector.x, thisVector.y, thisDir / PI * 180);
		//printf("\nAngle between vert %d (angle = %f) and %d (angle = %f) = %f", i, thisDir /PI *180, (i + 1) % 4, nextDir / PI * 180, thisAngleChange / PI * 180);
		//totalAngle = (theta < 0 ? totalAngle - thisAngle : totalAngle + thisAngle);
	}
	//printf("Angle = %f", totalAngle);
	printf("\nThe singularity at %f, %f is a ", v->x, v->y);
	if (totalAngle > -2*PI - .0001 && totalAngle < -2*PI + .0001) {
		printf("saddle point.\n\n");
		v->G = 1;
	}
	else {
		printf("not a saddle point.\n\n");
		v->B = 1;
	}

}


void drawSingularities() {
	for (auto v : Singularities) {
		drawDot(v->x, v->y, 0, .2, v->R, v->G, v->B);
	}
}



Vertex* getVertexAt(float xpos, float ypos) {
	
	if (xpos < min_x || xpos > max_x || ypos < min_y || ypos > max_y)
		return nullptr;

	Quad* thisFace;
	bool found = false;
	float minx , miny, maxx, maxy;
	int i;
	for (i = 0; i < poly->nquads; i++) {
		minx = miny = FLT_MAX;
		maxx = maxy = FLT_MIN;
		
		//check all the points for x and y to see if the point lies within this face
		for (int j = 0; j < 4; j++) {
			Vertex* thisVert = poly->qlist[i]->verts[j];
			if (thisVert->x < minx)
				minx = thisVert->x;
			if (thisVert->x > maxx)
				maxx = thisVert->x;
			if (thisVert->y < miny)
				miny = thisVert->y;
			if (thisVert->y > maxy)
				maxy = thisVert->y;
		}

		if (xpos >= minx && xpos <= maxx && ypos >= miny && ypos <= maxy) {
			//found correct face that contains this point
			thisFace = poly->qlist[i];
			found = true;
			goto next;
		}

	}
	next:
	if (found) {
		//now interpolate on vx, vy, and vz(?)
		/*printf("test point at %f , %f", xpos, ypos);
		printf("\nFound point in face %d:\n", i);
		printf("   vert 0: x,y: %f, %f , vx,vy: %f, %f \n", thisFace->verts[0]->x, thisFace->verts[0]->y, thisFace->verts[0]->vx, thisFace->verts[0]->vy);
		printf("   vert 1: x,y: %f, %f , vx,vy: %f, %f\n", thisFace->verts[1]->x, thisFace->verts[1]->y, thisFace->verts[1]->vx, thisFace->verts[1]->vy);
		printf("   vert 2: x,y: %f, %f , vx,vy: %f, %f\n", thisFace->verts[2]->x, thisFace->verts[2]->y, thisFace->verts[2]->vx, thisFace->verts[2]->vy);
		printf("   vert 3: x,y: %f, %f , vx,vy: %f, %f\n", thisFace->verts[3]->x, thisFace->verts[3]->y, thisFace->verts[3]->vx, thisFace->verts[3]->vy);
		//not sure if important, but vert 0 is max in x and y, 2 is min

		/*

			1 - 0
			|   |
			2 - 3

		*/



		//float vx = ()
		//interpolate in x direction first:
		float xLower = (maxx - xpos) / (maxx - minx) * thisFace->verts[2]->vx + (xpos - minx) / (maxx - minx) * thisFace->verts[3]->vx;
		float xUpper = (maxx - xpos) / (maxx - minx) * thisFace->verts[1]->vx + (xpos - minx) / (maxx - minx) * thisFace->verts[0]->vx;

		float yLower = (maxx - xpos) / (maxx - minx) * thisFace->verts[2]->vy + (xpos - minx) / (maxx - minx) * thisFace->verts[3]->vy;
		float yUpper = (maxx - xpos) / (maxx - minx) * thisFace->verts[1]->vy + (xpos - minx) / (maxx - minx) * thisFace->verts[0]->vy;

		//now interpolate in the y direction:


		float vxFinal = (maxy - ypos) / (maxy - miny) * xLower + (ypos - miny) / (maxy - miny) * xUpper;
		float vyFinal = (maxy - ypos) / (maxy - miny) * yLower + (ypos - miny) / (maxy - miny) * yUpper;

		//printf("\nestimation of the deltas at %f , %f\n", xpos, ypos);
		//printf("vx: %f, vy: %f", vxFinal, vyFinal);

		Vertex* v;
		v = new Vertex((double)xpos, (double)ypos, 0);
		v->vx = vxFinal;
		v->vy = vyFinal;
		return v;
	}
	return nullptr;
}





/******************************************************************************
Callback function for glut window reshaped
******************************************************************************/

void reshape(int width, int height) {

	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	/*Update pixels buffer*/
	free(pixels);
	pixels = (unsigned char *)malloc(sizeof(unsigned char)*win_width*win_height * 3);
	memset(pixels, 255, sizeof(unsigned char)*win_width*win_height * 3);
}


/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 5.0, 5.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected vert (quad) id = %d\n", poly->selected_vertex);
				
				foundSingleStreamline = false;

				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/*Display IBFV*/
void makePatterns(void)
{
	pixels = (unsigned char *)malloc(sizeof(unsigned char)*win_width*win_height * 3);
	memset(pixels, 255, sizeof(unsigned char)*win_width*win_height * 3);

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k, t;

	for (i = 0; i < 256; i++) 
		lut[i] = i < 127 ? 0 : 255;

	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) 
			phase[i][j] = rand() % 256;

	for (k = 0; k < Npat; k++) {
		t = k * 256 / Npat;
		for (i = 0; i < NPN; i++)
			for (j = 0; j < NPN; j++) {
				pat[i][j][0] =
					pat[i][j][1] =
					pat[i][j][2] = lut[(t + phase[i][j]) % 255];
				pat[i][j][3] = int(0.12 * 255);
			}
		glNewList(k + 1, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
		glEndList();
	}

}

void displayIBFV(void)
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_DEPTH_TEST);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/*draw the model with using the pixels, using vector field to advert the texture coordinates*/
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix1[16], projection_matrix1[16];
	int viewport1[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix1);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix1);
	glGetIntegerv(GL_VIEWPORT, viewport1);

	for (int i = 0; i < poly->nquads; i++) { //go through all the quads

		Quad *temp_q = poly->qlist[i];

		glBegin(GL_QUADS);

		for (int j = 0; j < 4; j++) {
			Vertex *temp_v = temp_q->verts[j];

			double x = temp_v->x;
			double y = temp_v->y;

			double tx, ty, dummy;

			gluProject((GLdouble)temp_v->x, (GLdouble)temp_v->y, (GLdouble)temp_v->z,
				modelview_matrix1, projection_matrix1, viewport1, &tx, &ty, &dummy);

			tx = tx / win_width;
			ty = ty / win_height;

			icVector2 dp = icVector2(temp_v->vx, temp_v->vy);
			normalize(dp);

			double dx = dp.x;
			double dy = dp.y;

			double r = dx * dx + dy * dy;
			if (r > dmax*dmax) {
				r = sqrt(r);
				dx *= dmax / r;
				dy *= dmax / r;
			}

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}

	iframe = iframe + 1;

	glEnable(GL_BLEND);

	/*blend the drawing with another noise image*/
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();


	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(iframe % Npat + 1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);


	/*draw the model with using pixels, note the tx and ty do not take the vector on points*/
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->nquads; i++) { //go through all the quads
		Quad *temp_q = poly->qlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++) {
			Vertex *temp_v = temp_q->verts[j];
			double x = temp_v->x;
			double y = temp_v->y;
			double tx, ty, dummy;
			gluProject((GLdouble)temp_v->x, (GLdouble)temp_v->y, (GLdouble)temp_v->z,
				modelview_matrix1, projection_matrix1, viewport1, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_BLEND);
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	CHECK_GL_ERROR();

	set_scene(GL_RENDER, poly);
	CHECK_GL_ERROR();

	/*display the mesh*/
	display_polyhedron(poly);
	CHECK_GL_ERROR();

	/*display selected elements*/
	display_selected_vertex(poly);
	CHECK_GL_ERROR();

	display_selected_quad(poly);
	CHECK_GL_ERROR();

	if(streamlines)
		drawPolyline(Streamlines, 1.5, 1, 0, 0);

	if(foundSingleStreamline)
		drawPolyline(SingleStreamline, 1.5, 0, 0, 1);

	if (singularities) 
		drawSingularities();
	

	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

/*global variable to save polylines*/
PolyLine pentagon;

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;
	case '4':
		display_mode = 4;
		{
			//examples for dot drawing and polyline drawing

			//create a polylines of a pentagon
			//clear current polylines
			pentagon.clear();
			//there are five vertices of a pentagon
			//the angle of each edge is 2pi/5.0
			double da = 2.0 * PI / 5.0;
			for (int i = 0; i < 5; i++) {
				double angle = i * da;
				double cx = cos(angle);
				double cy = sin(angle);

				double n_angle = (i + 1) % 5 * da;
				double nx = cos(n_angle);
				double ny = sin(n_angle);

				LineSegment line(cx, cy, 0, nx, ny, 0);
				pentagon.push_back(line);
			}

		}
		glutPostRedisplay();
		break;

	case '5':
		display_mode = 5;
		//show the IBFV of the field
		break;

	case 'r':
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;

	case '6':
		if (streamlines == false) {
			streamlines = true;
			if (foundStreamlines == false)
				drawAllStreamlines();
			foundStreamlines = true;
		}
		else {
			streamlines = false;
		}
		break;
	case '7':
		if (singularities == false) {
			singularities = true;
			if(foundSingularities == false)
				findSingularities();
			foundSingularities = true;
		}
		else {
			singularities = false;
		}
		break;
	}
}


/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/


void display_polyhedron(Polyhedron* poly)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	CHECK_GL_ERROR();

	switch (display_mode) {
	case 1:
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 1.0, 1.0, 0.0, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		CHECK_GL_ERROR();
	}
	break;

	case 2:
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();

		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
		break;

	case 4:
	{
		//draw a dot at position (0.2, 0.3, 0.4) 
		//with radius 0.1 in color blue(0.0, 0.0, 1.0)
		drawDot(0.2, 0.3, 0.4, 0.1, 0.0, 0.0, 1.0);

		//draw a dot at position of vlist[110]
		//with radius 0.2 in color magenta (1.0, 0.0, 1.0)
		Vertex *v = poly->vlist[110];
		drawDot(v->x, v->y, v->z, 0.2, 1.0, 0.0, 1.0);

		//draw line segment start at vlist[110] and end at (vlist[135]->x, vlist[135]->y, 4)
		//with color (0.02, 0.1, 0.02) and width 1
		LineSegment line(poly->vlist[110]->x, poly->vlist[110]->y, poly->vlist[110]->z,
			poly->vlist[135]->x, poly->vlist[135]->y, 4);
		drawLineSegment(line, 1.0, 0.0, 1.0, 0.0);

		//draw a polyline of pentagon with color orange(1.0, 0.5, 0.0) and width 2
		drawPolyline(pentagon, 2.0, 1.0, 0.5, 0.0);

		//display the mesh with color cyan (0.0, 1.0, 1.0)
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(0.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 5:
		displayIBFV();
		break;
	}
}
