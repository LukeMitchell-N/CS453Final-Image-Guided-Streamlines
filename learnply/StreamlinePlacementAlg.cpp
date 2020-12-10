/*
Implementation of Turk and Banks' image-based streamline placement algorithm

This algorithm takes an initial set of streamlines of a vector field and iteratively improves upon them
The algorithm uses a filtered image of the streamlines to create an energy value for the set, which is used to compare to other potential sets to find the best placed streamlines
*/

#include "polyline.h"

#define resolution 128
float image[resolution][resolution];
float targetEnergy = .5;
float radiusOfInfluence = 1.;
double minXY = -10., maxXY = 10.;			//taking for granted that all the ply files have the same bounds

//forward declarations of functions
void optimizeStreamlines(StreamlineSet*);
void streamlinesToImage(StreamlineSet*);
void streamlineToImage(PolyLine*);

void visualizeImage();

//iteratively alter the streamline set
void optimizeStreamlines(StreamlineSet* set) {

	int d = set->size();
	streamlinesToImage(set);

	visualizeImage();

}


void streamlinesToImage(StreamlineSet* set) {

	for (int i = 0; i < set->size(); i++) {
		streamlineToImage(set->at(i).p);
	}

}

void streamlineToImage(PolyLine* pl) {

	int gridRadius = radiusOfInfluence / (maxXY - minXY) * resolution;
	for (int i = 0; i < pl->size(); i++) {
		int gridx = (pl->at(i).start.x - minXY) / (maxXY - minXY) * resolution;
		int gridy = (int)(pl->at(i).start.y - minXY) / (maxXY - minXY) * resolution;

		//because we know the radius of influence for each point, there is no need to calculate the influence outside of the range
		for (int j = -gridRadius; j < gridRadius; j++) {
			for (int k = -gridRadius; k < gridRadius; k++) {
				if (gridy + j > 0 && gridy + j < resolution && gridx + k > 0 && gridx + k < resolution) {

					//using Turk and Banks' Kernel function
					float r = sqrt(j * j + k * k) / gridRadius;
					if (r < 1)
						image[gridy + j][gridx + k] += 2 * (r * r * r) - 3 * (r * r) + 1;
				}
			}
		}


	}

}

//very crude visualization of the filtered image
void visualizeImage() {
	for (int i = 0; i < resolution; i++) {
		for (int j = 0; j < resolution; j++){
			if (image[i][j] == 0)
				printf(".");
			else if (image[i][j] >= 5)
				printf("#");
			else
				printf("\"");
		}
		printf("\n");
	}
}


