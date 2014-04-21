#ifndef __OBJECT_H__
#define __OBJECT_H__



#define OBJ_CUBE						10
#define OBJ_TRUNCATED_ICOSAHEDRON		11
#define OBJ_TETRAHEDRON					12
#define OBJ_DODECAHEDRON				13

//#include <vector>
#include "../defines.hpp"
#include <cstring>

struct Edge;

struct Face {
	Edge *from;		// Första edgen
	Vec Norm;		// Face Normal
};

struct Edge {
	Vec *fr;
	Vec *to;
	Edge *next;
	Edge *prev;
	Edge *oppo;
	Face *face;

	void set(Vec* _V, Edge *_E, Face *_F, int _fr, int _to, int _next, int _prev, int _oppo, int _face);
};



	// Object med facenormals
class ObjectFN {
protected:

	int numV;
	int numE;
	int numF;

	Vec *V;
	Edge *E;
	Face *F;	// 



public:
	ObjectFN();
	~ObjectFN();
	ObjectFN(int _vert, int _edge, int _face);
	ObjectFN(std::string);

	void print();

	Face* getFaces(int &numOfFaces) const;
};



class ObjCubeFN: public ObjectFN{
public:
	ObjCubeFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

class ObjTetrahedronFN: public ObjectFN{
public:
	ObjTetrahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};

//OBJ_DODECAHEDRON
class ObjDodecahedronFN: public ObjectFN {
public:
	ObjDodecahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori);
};


class World {
private:
	static std::list<ObjectFN*> Objs;
	
public:
	static bool addObjectFN(int objType, const Vec &Pos, const Vec &Siz, const Mat &Ori);
	static bool removeAllObjects();
	static ObjectFN* getAnObject();
};



#endif 
