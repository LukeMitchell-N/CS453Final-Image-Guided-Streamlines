# CS453Final-Image-Guided-Streamlines
Implementation of Turk and Banks' 1994 Image Guided Streamline Placement algorithm
http://cs.brown.edu/courses/cs237/2000/1999/turk.pdf

~~~~ Design ~~~~

The algorithm by Turk and Banks uses a Low-Pass filtered image of the streamline set for comparison between streamline sets
  This can be stored as a 2d array of values from 0 to 1 or any range
  Store the dimensions as globals, in case we want to create a high resolution image for possible output
  We could use two 2d arrays, a main an a test, we set changes in test, and if we prefer the changes, we set the main to point to test vice versa
  
Streamlines
  Streamlines need to be drawn with determinable length or number of segments
  Each streamline (class or struct) will need to include ties to information such as:
    The locations of the endpoints of its segments (this will aid in calculations)
    Some method of storing (vector seems promising) which pixels in the blurred image this point affects and by how much
  
 The effect of each streamline point on the blurred image will need to be it's own function
  The energy produced by a point is found by a circularly symetric function that falls off further from the point
  The energy a point applies to a given pixel is = 2r3 - 3r2 +1. This  can be found on the bottom of page 3 of the paper 
  
