#include "objSimple.hpp"

ObjSimpleFN::ObjSimpleFN(const Vec &Pos, const Vec &Siz, const Mat &Ori) :ObjectFN(3, 6, 2)
{
	//	cout << "Ska skapa en Kub" << endl;
	V[0].X = Vec(1,			0,			0);		V[0].from = &E[0];
	V[1].X = Vec(-.5, 		sqrt(.75),	0);		V[1].from = &E[1];
	V[2].X = Vec(-.5,		-sqrt(.75),	0);		V[2].from = &E[2];

	for (int i=0; i<3; i++)
		V[i].X = Vec(V[i].X.x*Siz.x, V[i].X.y*Siz.y, V[i].X.z*Siz.z) * Ori + Pos;

	E[0].set(V, E, F, 0, 1, 1, 2, 3, 0);
	E[1].set(V, E, F, 1, 2, 2, 0, 4, 0);
	E[2].set(V, E, F, 2, 0, 0, 1, 5, 0);

	E[3].set(V, E, F, 1, 0, 5, 4, 0, 1);
	E[4].set(V, E, F, 2, 1, 3, 5, 1, 1);
	E[5].set(V, E, F, 0, 2, 4, 3, 2, 1);


	F[0].from = &E[0];		F[0].Norm = Vec(0, 0, 1) * Ori;
	F[1].from = &E[3];		F[1].Norm = Vec(0, 0, -1) * Ori;

	consistsOfOnlyTriangles = true;
}