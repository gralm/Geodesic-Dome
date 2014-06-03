#ifndef __OBJECT_H__
#define __OBJECT_H__



#define OBJ_CUBE						10
#define OBJ_TRUNCATED_ICOSAHEDRON		11
#define OBJ_TETRAHEDRON					12
#define OBJ_DODECAHEDRON				13

// https://www8.cs.umu.se/kurser/5DV009/HT08/handouts/HO%20E%20-%20Subdivision.pdf
#include <vector>
#include "../defines.hpp"
#include <cstring>

#define _HI(v)		(((v)+1) % 5)
#define _LO(v)		((v)? (v-1): (4))



struct Edge;

struct Vertex {
	Edge *from;
	Vec X;
	Vec Norm;

	void print(const Vertex *zero);
};

struct Face {
	Edge *from;		// Första edgen
	Vec Norm;		// Face Normal

	void print(const Face *zero);
	void update();
};

struct Edge {
	Vertex *fr;
	Vertex *to;
	Edge *next;
	Edge *prev;
	Edge *oppo;
	Face *face;

	void set(Vertex* _V, Edge *_E, Face *_F, int _fr, int _to, int _next, int _prev, int _oppo, int _face);
	void print(const Edge *zero);
};




	// Object med facenormals
class ObjectFN {
protected:

	bool consistsOfOnlyTriangles;

	int numV;
	int numE;
	int numF;

	//Vec *V;
	Vertex *V;
	Edge *E;
	Face *F;	// 



public:
	ObjectFN();
	~ObjectFN();
	ObjectFN(int _vert, int _edge, int _face);
	ObjectFN(std::string);

	bool subdivide1();
	bool subdivide2();

	bool updateConsistsOfOnlyTriangles();

		// http://en.wikipedia.org/wiki/Dual_polyhedron
	bool makeDual();
	//ObjectFN *createDual();
	ObjectFN *truncate(TYP val);	// truncated = 0.5, rectified = 1.0;
	ObjectFN *rectify();
	ObjectFN *truncate();


	TYP normalizeRadius();
	bool test() const;

	void print();

	Face* getFaces(int &numOfFaces) const;
};


class World {
private:
	static std::list<ObjectFN*> Objs;
	
public:
	static ObjectFN* addObjectFN(int objType, const Vec &Pos, const Vec &Siz, const Mat &Ori);
	static bool addObjectFN(ObjectFN *pObj);
	static bool removeAllObjects();
	static ObjectFN* getAnObject();
	static std::list<ObjectFN*>* getObjectListPointer();
};


#endif 
