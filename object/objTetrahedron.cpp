#include "objTetrahedron.hpp"

ObjTetrahedronFN::ObjTetrahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori): ObjectFN(4, 12, 4)
{
	TYP s2 = static_cast<TYP>(sqrt(2.0));
	TYP s3 = static_cast<TYP>(sqrt(3.0));

	V[0].X = Vec(1/s3, 		0, 		-0.5/(s3*s2));
	V[1].X = Vec(-0.5/s3, 	0.5, 	-0.5/(s3*s2));
	V[2].X = Vec(-0.5/s3, 	-0.5, 	-0.5/(s3*s2));
	V[3].X = Vec(0, 		0, 		(s3*s2)*.25);
	
	V[0].from = &E[0];
	V[1].from = &E[1];
	V[2].from = &E[4];
	V[3].from = &E[2];

	for (int i=0; i<4; i++)
		V[i].X = Vec(V[i].X.x*Siz.x, V[i].X.y*Siz.y, V[i].X.z*Siz.z) * Ori + Pos;


	E[0].set(V, E, F, 0, 1, 1, 2,11, 0);
	E[1].set(V, E, F, 1, 3, 2, 0, 8, 0);
	E[2].set(V, E, F, 3, 0, 0, 1, 5, 0);
	E[3].set(V, E, F, 3, 2, 4, 5, 7, 1);
	E[4].set(V, E, F, 2, 0, 5, 3, 9, 1);
	E[5].set(V, E, F, 0, 3, 3, 4, 2, 1);
	E[6].set(V, E, F, 1, 2, 7, 8,10, 2);
	E[7].set(V, E, F, 2, 3, 8, 6, 3, 2);
	E[8].set(V, E, F, 3, 1, 6, 7, 1, 2);
	E[9].set(V, E, F, 0, 2,10,11, 4, 3);
	E[10].set(V,E, F, 2, 1,11, 9, 6, 3);
	E[11].set(V,E, F, 1, 0, 9,10, 0, 3);

	F[0].from = &E[0];		F[0].Norm = (V[2].X * (-2.0*s2/s3)) * Ori;
	F[1].from = &E[3];		F[1].Norm = (V[1].X * (-2.0*s2/s3)) * Ori;
	F[2].from = &E[6];		F[2].Norm = (V[0].X * (-2.0*s2/s3)) * Ori;
	F[3].from = &E[9];		F[3].Norm = (V[3].X * (-2.0*s2/s3)) * Ori;


	//subdivide1();
}

