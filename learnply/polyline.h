#pragma once
#include "icVector.H"
#include <vector>

class LineSegment
{
public:
	LineSegment(icVector3 s, icVector3 e) {
		start = s;
		end = e;
		len = length(end - start);
	};
	LineSegment(double sx, double sy, double sz, double ex, double ey, double ez) {
		start = icVector3(sx, sy, sz);
		end = icVector3(ex, ey, ez);
		len = length(end - start);
	};
	LineSegment() {
		start.x = 0;
		end.y = 0;
		len = 0;
	};
	
	icVector3 start, end;  // all in local coordinate systems
	double len;
};

typedef std::vector<LineSegment> PolyLine;
