# CS453Final-Image-Guided-Streamlines
Implementation of Turk and Banks' 1994 Image Guided Streamline Placement algorithm
http://cs.brown.edu/courses/cs237/2000/1999/turk.pdf


Design:
~~~ 

The algorithm by Turk and Banks uses a Low-Pass filtered image of the streamline set for comparison between streamline sets
  This can be stored as a 2d array of values from 0 to 1 or any range
  Store the dimensions as globals, in case we want to create a high resolution image for possible output
  We could use two 2d arrays, a main an a test, we set changes in test, and if we prefer the changes, we set the main to point to test vice versa
    (thought more about this, we might not need this if we can just check the energy of a specific streamline)
  
Streamlines
  Streamlines need to be drawn with determinable length or number of segments
  Each streamline (class or struct) will need to include ties to information such as:
    The locations of the endpoints of its segments (this will aid in calculations)
    Some method of storing (vector seems promising) which pixels in the blurred image this point affects and by how much
  
The effect of each streamline point on the blurred image will need to be it's own function
  The energy produced by a point is found by a circularly symetric function that falls off further from the point
  The energy a point applies to a given pixel is = 2r3 - 3r2 +1. This  can be found on the bottom of page 3 of the paper 
  
Operations and picking streamlines
  The paper outlines what they call an "oracle", which maintains a priority queue for which streamlines contribute the most to the energy
  We can determine how much effort we want to put into the oracle (the paper describes ways in which it reads the image data and can suggest new inserts, and it recommends specific changes for streamlines, but we could do purely random changes if we wanted)
  We could have the algorithm pick a random streamline and a random operation (delete, shorten, lengthen, move origin)
~~~

Here are some of my thoughts for tasks:

- Make/ensure that streamlines can be drawn for a set number of steps/length forward and back, from an origin point
- Make the streamlines draw in some way (I can't remember what my code did for project 3. We should make the seeds a grid)
- Make a streamline class/struct and give it the properties listed above (needs to be done early and with some thought)
- Make global variables for the image (2d Array of floats)

Functions:
- Determine energy of full image
- Determine energy from a single streamline
- Draw the entire set of streamlines to the screen
- Iterative function that determines a new random change and calls the right function
- Move seed point function
- Shorten streamline function
- Lengthen streamline function
- Delete streamline function
- Add streamline function ?
- Merge streamlines functions ?
