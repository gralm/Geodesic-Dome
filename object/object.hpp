#ifndef __OBJECT_H__
#define __OBJECT_H__

#define OBJ_DIR_CW						0
#define OBJ_DIR_CCW						1


#define OBJ_CUBE						10
#define OBJ_TRUNCATED_ICOSAHEDRON		11
#define OBJ_TETRAHEDRON					12
#define OBJ_DODECAHEDRON				13
#define OBJ_TRUNCATED_TETRAHEDRON		14
#define OBJ_TRUNCATED_OCTAHEDRON		15
#define OBJ_ICOSAHEDRON 				16
#define OBJ_OCTAHEDRON 					17
#define OBJ_TRUNCATED_CUBE				18
#define OBJ_SIMPLE						19
#define OBJ_ICOSIDODECAHEDRON			20
#define OBJ_CUBOCTAHEDRON				21
#define OBJ_RHOMBICOSIDODECAHEDRON		22
#define OBJ_SNUBDODECAHEDRON			23
#define OBJ_SNUBDODECAHEDRON_CW			23
#define OBJ_SNUBDODECAHEDRON_CCW		24
#define OBJ_RHOMBICUBOCTAHEDRON			25
#define OBJ_SNUBCUBE					26
#define OBJ_SNUBCUBE_CW					26
#define OBJ_SNUBCUBE_CCW				27
#define OBJ_TRUNCATED_RHOMBIC_DODECAHEDRON	28
#define OBJ_CHAMFERED_CUBE				28
#define OBJ_TRUNCATED_RHOMBIC_TRIACONTAHEDRON	29
#define OBJ_CHAMFERED_DODECAHEDRON		OBJ_TRUNCATED_RHOMBIC_TRIACONTAHEDRON

#define OBJ_GOLDBERG_0_1				OBJ_DODECAHEDRON
#define OBJ_GOLDBERG_1_1				OBJ_TRUNCATED_TETRAHEDRON
#define OBJ_GOLDBERG_1_2				1012
#define OBJ_GOLDBERG_0_2				OBJ_CHAMFERED_DODECAHEDRON




#define OBJ_TYPE_CONCAVE				1
#define OBJ_TYPE_CONVEX					2
#define OBJ_TYPE_IRREGULAR				4
#define OBJ_TYPE_REGULAR 				8
#define OBJ_TYPE_SPHERICAL				0x12			// sfär = 0x10 + att den är konvex 

// https://www8.cs.umu.se/kurser/5DV009/HT08/handouts/HO%20E%20-%20Subdivision.pdf
#include <vector>

#include <list>
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
	int countEdges();
};

struct Face {
	Edge *from;		// Första edgen
	Vec Norm;		// Face Normal

	void print(const Face *zero);
	void update();
	Vec getCenter();
	int countEdges();
	void moveUp(TYP height);
	TYP maxSinErr(int &N);		// , N = num of edges
	void rotate(TYP rad);
};

struct Edge {
	Vertex *fr;
	Vertex *to;
	Edge *next;
	Edge *prev;
	Edge *oppo;
	Face *face;

	void set(Vertex* _V, Edge *_E, Face *_F, int _fr, int _to, int _next, int _prev, int _oppo, int _face);
	void print(Vertex *V0, Edge *E0, Face *F0);
	TYP length();
};



	// Object med facenormals
class ObjectFN {
protected:

	bool consistsOfOnlyTriangles;

	int numV;
	int numE;
	int numF;

	Vertex *V;
	Edge *E;
	Face *F;	// 

	bool CopyVEF(Vertex *nyV, Edge *nyE, Face *nyF);	// Kopiera alla som inte är 0
	bool EmptyVEF(Vertex *nyV, Edge *nyE, Face *nyF);
	static int truncatedEdgeNum(int _N, int _r, int _n, int _p);
	static Vertex *truncatedVertexNum(Vertex *_Start, int _N, int _r, int _n, int _p);
	static Edge *truncatedEdgeNum(Edge *_Start, int _N, int _r, int _n, int _p);
	static void print(const Vertex *V_, const Edge *E_, const Face *F_, int numV_, int numE_, int numF_);

public:
	ObjectFN();
	~ObjectFN();
	ObjectFN(int _vert, int _edge, int _face);
	ObjectFN(std::string);

	void edgeCompare(TYP&, TYP&);

	ObjectFN *greenHousify(TYP b, TYP h);

	bool subdivide1();		// makes pyramids of all surfaces
	bool subdivide1(int n, TYP height);	// makes pyramids of surfaces with n edges
	
	bool subdivide2();		// subdivides triangles to four trianlges
	bool subdivide2(int n);	// divides every edge n times. subdivide2() = subdivide(2)

	bool updateConsistsOfOnlyTriangles();

		// http://en.wikipedia.org/wiki/Dual_polyhedron
	bool makeDual();
	bool truncate(TYP val);		// 0 < val < 1,		0 < truncated < 1, rectified = 1.0;

	bool truncate2(TYP val);	// 0 < val < 1,		truncated = 0.5, rectified = 1.0;
	bool rectify();				// truncate(val = 1.0)
	//bool snub(int n);

		//nya egdes får längden: val * sqrt(2 - 2*NA*NB),  där NA och NB är två grannsidors normaler.
	bool expand(TYP val);		// val är en radiella förändringsfaktorn, val > 1.
	bool chamfer(TYP val);		// skapar hexagoner mellan edges, alla hörnpunkt i polyedern måste vara ansluten till 3 edges.
	bool rotatePolygons(TYP angle, int N);	// Rotate faces with N vertices by angle
	bool snub(TYP h, TYP t, int N);
	bool splitBrokenTetragons();
	bool setPolygonHeight(TYP h, int N);	// sätt höjden, h, på vertiecar från centrum på polygoner med N edges

	Vec getCenter();

		// Ange noll om du inte vill göra
	bool transform(const Vec *Pos, const Vec *Siz, const Mat *Ori);	

	TYP normalizeRadius();
	bool normalizeNormals();
	bool test(unsigned int shapeType_) const;

	void print();
	Face* getFaces(int &numOfFaces) const;


	void tabort();

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
