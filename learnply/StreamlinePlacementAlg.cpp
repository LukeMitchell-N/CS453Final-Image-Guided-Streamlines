/*
Implementation of Turk and Banks' image-based streamline placement algorithm

This algorithm takes an initial set of streamlines of a vector field and iteratively improves upon them
The algorithm uses a filtered image of the streamlines to create an energy value for the set, which is used to compare to other potential sets to find the best placed streamlines
*/
