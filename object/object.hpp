#ifndef __OBJECT_H__
#define __OBJECT_H__

//#include <vector>
#include "../defines.hpp"
#include <cstring>

struct Edge;

struct Face {
	Edge *from;		// Första edgen
	Vec *Norm;		// Face Normal
};

struct Edge {
	Vec *fr;
	Vec *to;
	Edge *next;
	Edge *prev;
	Edge *oppo;
	Face *face;
};

	// Object med facenormals
class ObjectFN {
protected:

	//std::vector<Vec> V;
	//std::vector<Edge> E;
	Vec *V;
	Edge *E;
	Face *F;	// 


public:
	ObjectFN();
	ObjectFN(std::string);
};


#endif 