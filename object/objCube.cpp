#include "objCube.hpp"

ObjCubeFN::ObjCubeFN(const Vec &Pos, const Vec &Siz, const Mat &Ori) :ObjectFN(8, 24, 6)
{
	//	cout << "Ska skapa en Kub" << endl;
	V[0].X = Vec(-.5,		.5, 	-.5);
	V[1].X = Vec(.5, 		.5,	 	-.5);
	V[2].X = Vec(-.5,		.5, 	.5);
	V[3].X = Vec(.5, 		.5, 	.5);
	V[4].X = Vec(-.5,		-.5,	-.5);
	V[5].X = Vec(.5, 		-.5,	-.5);
	V[6].X = Vec(-.5,		-.5,	.5);
	V[7].X = Vec(.5, 		-.5,	.5);


	for (int i=0; i<8; i++)
		V[i].X = Vec(V[i].X.x*Siz.x, V[i].X.y*Siz.y, V[i].X.z*Siz.z) * Ori + Pos;

	E[0].set(V, E, F, 1, 0, 1, 3, 4, 0);
	E[1].set(V, E, F, 0, 2, 2, 0, 11, 0);
	E[2].set(V, E, F, 2, 3, 3, 1, 23, 0);
	E[3].set(V, E, F, 3, 1, 0, 2, 16, 0);

	E[4].set(V, E, F, 0, 1, 5, 7, 0, 1);
	E[5].set(V, E, F, 1, 5, 6, 4, 19, 1);
	E[6].set(V, E, F, 5, 4, 7, 5, 12, 1);
	E[7].set(V, E, F, 4, 0, 4, 6, 8, 1);

	E[8].set(V, E, F, 0, 4, 9, 11, 7, 2);
	E[9].set(V, E, F, 4, 6, 10, 8, 15, 2);
	E[10].set(V, E, F, 6, 2, 11, 9, 20, 2);
	E[11].set(V, E, F, 2, 0, 8, 10, 1, 2);

	E[12].set(V, E, F, 4, 5, 13, 15, 6, 3);
	E[13].set(V, E, F, 5, 7, 14, 12, 18, 3);
	E[14].set(V, E, F, 7, 6, 15, 13, 21, 3);
	E[15].set(V, E, F, 6, 4, 12, 14, 9, 3);

	E[16].set(V, E, F, 1, 3, 17, 19, 3, 4);
	E[17].set(V, E, F, 3, 7, 18, 16, 22, 4);
	E[18].set(V, E, F, 7, 5, 19, 17, 13, 4);
	E[19].set(V, E, F, 5, 1, 16, 18, 5, 4);

	E[20].set(V, E, F, 2, 6, 21, 23, 10, 5);
	E[21].set(V, E, F, 6, 7, 22, 20, 14, 5);
	E[22].set(V, E, F, 7, 3, 23, 21, 17, 5);
	E[23].set(V, E, F, 3, 2, 20, 22, 2, 5);


	F[0].from = &E[0];		F[0].Norm = Vec(0, 1, 0) * Ori;
	F[1].from = &E[4];		F[1].Norm = Vec(0, 0, -1) * Ori;
	F[2].from = &E[8];		F[2].Norm = Vec(-1, 0, 0) * Ori;
	F[3].from = &E[12];		F[3].Norm = Vec(0, -1, 0) * Ori;
	F[4].from = &E[16];		F[4].Norm = Vec(1, 0, 0) * Ori;
	F[5].from = &E[20];		F[5].Norm = Vec(0, 0, 1) * Ori;
	/*cout << "Skapade en kub" << endl;

	for (int i=0; i<8; i++)
		cout << V[i].X << endl;
	cout << "antalet killar \n";
	print();
	subdivide1();*/
}
