/*
Implementation of Turk and Banks' image-based streamline placement algorithm

Authors:	Luke Mitchell	- mitchluk
			Mason Hall		- hallmas

This algorithm takes an initial set of streamlines of a vector field and iteratively improves upon them
The algorithm uses a filtered image of the streamlines to create an energy value for the set, which is used to compare to other potential sets to find the best placed streamlines
*/

#include "polyline.h"
#include <time.h>

//image variables
#define			resolution 128
float			buffer1[resolution][resolution], buffer2[resolution][resolution];
float			(*image)[resolution][resolution] = &buffer1;
float			(*test)[resolution][resolution] = &buffer2;
float			(*swap)[resolution][resolution];

//image energy variables
#define			sufficientNumFails 200
#define			upperEnergyLimit 50
float			targetEnergy = 5;
float			radiusOfInfluence = 1;
double			minXY = -10., maxXY = 10.;			//taking for granted that all the ply files have the same bounds
const int		subDivisions = 5;

//streamline variables
StreamlineSet	streamlines;
StreamlineSet	testLines;
PolyLine		SingleStreamline;
bool			streamlinesOn = false, singularities = false;
bool			foundStreamlines = false, foundSingleStreamline = false, foundSingularities = false;
int				streamLength = 50;

//operation variables
#define randMaxMovement 1


//forward declarations of functions
void optimizeStreamlines(StreamlineSet*);
void streamlinesToImage(StreamlineSet*, float[][resolution]);
void streamlineAlterImage(PolyLine*, float[][resolution], bool);
float getEnergyFromImage(float[][resolution]);
void copyImage(float[][resolution], float[][resolution]);
void visualizeImage(float[][resolution]);
bool areStreamlinesSatisfactory(int, float[][resolution]);

//operations
//note: all operations implicitly operate on the test image, which is compared after the operation is complete
void adjustRandomSeed(int);


//streamline functions
void placeStreamlinesGrid();
void drawLineRecursive(Vertex* v, PolyLine*, bool forward, int);
Vertex* getVertexAt(float xpos, float ypos);
Vertex* RKGetNextVertex(Vertex* initial, bool forward);




//iteratively alter the streamline set
void optimizeStreamlines(StreamlineSet* set) {


	int d = set->size();
	streamlinesToImage(set, *image);
		visualizeImage(*image);
		float e = getEnergyFromImage(*image);
		float normalizedEnergy = e / (resolution * resolution);
		//printf("This image has a normalized energy value of %f", normalizedEnergy);

	srand(time(NULL));
	int numFails = 0;
	float baseEnergy, testEnergy;
	while (areStreamlinesSatisfactory(numFails, *image) == false) {
		//randomly select an operation to perform on the streamlines
		//but for now move a random seed point
		int randStreamline = rand() % set->size();
		copyImage(*image, *test);
		adjustRandomSeed(randStreamline);
		
		baseEnergy = getEnergyFromImage(*image);
		testEnergy = getEnergyFromImage(*test);
		printf("Old line placement resulted in %f energy, new line placement resulted in %f energy\n", baseEnergy, testEnergy);
		
		//if this is an improvement, sway the buffers
		if (testEnergy < baseEnergy) {
			swap = image;
			image = test;
			test = swap;
		}
	}
}

//write all streamlines to image
void streamlinesToImage(StreamlineSet* set, float img[][resolution]) {

	for (int i = 0; i < set->size(); i++) {
		streamlineAlterImage(set->at(i).p, img, true);
	}

}

//function to alter an image based on the points of a streamline.
//Provides functionality to add influence to an image as well as take influence away
void streamlineAlterImage(PolyLine* pl, float img[][resolution], bool addInfluence) {

	int gridRadius = radiusOfInfluence / (maxXY - minXY) * resolution;
	for (int i = 0; i < pl->size(); i++) {
		//printf("streamline point %d is at point %f, %f \n", i, pl->at(i).start.x, pl->at(i).start.y);
		int gridx = (pl->at(i).start.x - minXY) / (maxXY - minXY) * resolution;
		int gridy = (int)(pl->at(i).start.y - minXY) / (maxXY - minXY) * resolution;

		//because we know the radius of influence for each point, there is no need to calculate the influence outside of the range
		for (int j = -gridRadius; j < gridRadius; j++) {
			for (int k = -gridRadius; k < gridRadius; k++) {
				//printf("calculating image effect of point %d\n", ++testCounter);
				if (gridy + j > 0 && gridy + j < resolution && gridx + k > 0 && gridx + k < resolution) {

					//using Turk and Banks' Kernel function
					float r = sqrt(j * j + k * k) / gridRadius;
					if (r < 1) {
						if (addInfluence)
							img[gridy + j][gridx + k] += 2 * (r * r * r) - 3 * (r * r) + 1;		//either add influence
						else
							img[gridy + j][gridx + k] -= 2 * (r * r * r) - 3 * (r * r) + 1;		//or remove it

					}
					//if (img[gridy + j][gridx + k] > 1)										//cap each value to 1
						//img[gridy + j][gridx + k] = 1;
				}
			}
		}


	}

}


//copy one image buffer to another
void copyImage(float src[][resolution], float dest[][resolution]) {
	for (int i = 0; i < resolution; i++) {
		for (int j = 0; j < resolution; j++) {
			dest[i][j] = src[i][j];
		}
	}
}

//very crude visualization of the filtered image
void visualizeImage(float img[][resolution]){
	for (int i = 0; i < resolution; i++) {
		for (int j = 0; j < resolution; j++){
			if (img[i][j] == 0)
				printf(".");
			else if (img[i][j] >= 30)
				printf("#");
			else
				printf("\"");
		}
		printf("\n");
	}
}


float getEnergyFromImage(float img[][resolution]) {
	float runningTotal = 0;
	
	for (int i = 0; i < resolution; i++) {
		for (int j = 0; j < resolution; j++) {
			runningTotal += (img[i][j] - targetEnergy) * (img[i][j] - targetEnergy);
		}
	}
	return runningTotal;
}


bool areStreamlinesSatisfactory(int numFails, float img[][resolution]) {
	if (numFails > sufficientNumFails || getEnergyFromImage(img) < targetEnergy * 1.2)
		return true;
	return false;
}




//operation functions

void adjustRandomSeed(int sl) {

	streamlineAlterImage(streamlines.at(sl).p,  *test,  false);			//first remove the streamline's influence on the image
	float newRandX = streamlines.at(sl).seed.x + (rand() % 100) / 100 * (randMaxMovement * 2) - randMaxMovement;
	float newRandY = streamlines.at(sl).seed.y + (rand() % 100) / 100 * (randMaxMovement * 2) - randMaxMovement;
	Vertex* v = getVertexAt(newRandX, newRandY);
	Streamline alteredLine;
	alteredLine.seed.x = newRandX;
	alteredLine.seed.y = newRandY;
	alteredLine.length = streamLength * 2;
	alteredLine.p = new PolyLine(streamLength * 2);
	drawLineRecursive(v, alteredLine.p, true, streamLength * 2);

	//finally add the influence of the new line to the test image and return
	streamlineAlterImage(alteredLine.p, *test, true);
	
}



//******************  Streamline placement functions *********************

void placeStreamlinesGrid() {

	float x, y;
	for (int i = 0; i < subDivisions; i++) {
		x = minXY + (1 + i) * (1. / ((float)subDivisions + 2)) * (maxXY - minXY);
		for (int j = 0; j < subDivisions; j++) {
			y = minXY + (1 + j) * (1. / ((float)subDivisions + 2)) * (maxXY - minXY);
			Vertex* v = getVertexAt(x, y);
			streamlines.at(i * subDivisions + j).length = streamLength * 2;
			streamlines.at(i * subDivisions + j).seed.x = v->x;
			streamlines.at(i * subDivisions + j).seed.y = v->y;
			streamlines.at(i * subDivisions + j).p = new PolyLine(streamLength *2);
			drawLineRecursive(v, streamlines.at(i * subDivisions + j).p,  true, streamLength);
			drawLineRecursive(v, streamlines.at(i * subDivisions + j).p, false, streamLength);
		}
	}
	foundStreamlines = true;

}

void drawLineRecursive(Vertex* v, PolyLine* p, bool forward, int depthRemaining) {

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

	//if this streamline is going in the backward direction, reverse the start and end
	if (forward == false) {
		newLine.start = { next_v->x, next_v->y, 0 };
		newLine.end = { v->x, v->y, 0 };
	}

	p->push_back(newLine);


	//continue the line from the next point
	drawLineRecursive(next_v, p, forward, --depthRemaining);

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
	double timestep = .2 * direction, timeStepProportion;		//I think this isn't enough to make the forward/back work right, revisit later
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




Vertex* getVertexAt(float xpos, float ypos) {

	if (xpos < minXY || xpos > maxXY || ypos < minXY || ypos > maxXY)
		return nullptr;

	Quad* thisFace;
	bool found = false;
	float minx, miny, maxx, maxy;
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

