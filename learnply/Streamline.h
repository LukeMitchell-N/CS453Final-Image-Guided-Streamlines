#include "Polyhedron.h"
#include <vector>


class Streamline {

	
public:
	icVector2 seed;
	int length;
	PolyLine *p;
	Streamline() {
		p = new PolyLine;
	}

private:


};

typedef std::vector<Streamline> StreamlineSet;