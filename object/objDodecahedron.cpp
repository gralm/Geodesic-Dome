#include "objDodecahedron.hpp"
	//OBJ_DODECAHEDRON
ObjDodecahedronFN::ObjDodecahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori): ObjectFN(20, 60, 12)
{
	TYP s5 = static_cast<TYP>(sqrt(5.0));
	TYP sin1 = sqrt((5.0+s5) / 8.);
	TYP sin2 = sqrt((5.0-s5) / 8.);

		// A
	V[0].X = Vec(1.0, 0.0, 0.0);
	V[1].X = Vec(0.25*(s5-1.0), sin1, 0.0);
	V[2].X = Vec(-0.25*(s5+1.0), sin2, 0.0);
	V[3].X = Vec(V[2].X.x, -V[2].X.y, 0.0);
	V[4].X = Vec(V[1].X.x, -V[1].X.y, 0.0);

		// B
	for (int i=0; i<5; i++)
		V[5+i].X = V[i].X*(0.5*(1+s5)) + Vec(0.0, 0.0, 0.25*(1.0-s5));
		//V[5+i] = V[i]*(0.25*(1+s5)) + Vec(0.0, 0.0, 0.25*(1.0-s5));

		// C
	for (int i=0; i<5; i++)
		V[10+i].X = -V[i<2? i+8: i+3].X;

		// A
	for (int i=0; i<5; i++)
		V[i].X.z -= 0.25*(3+s5);

		// D
	for (int i=0; i<5; i++)
		V[15+i].X = -V[(i<2)? i+3: i-2].X;


	E[0].set(V, E, F, 0, 4, 1, 4,25, 0);
	E[1].set(V, E, F, 4, 3, 2, 0,20, 0);
	E[2].set(V, E, F, 3, 2, 3, 1,15, 0);
	E[3].set(V, E, F, 2, 1, 4, 2,10, 0);
	E[4].set(V, E, F, 1, 0, 0, 3, 5, 0);



	E[5].set(V, E, F, 0, 1, 6, 9, 4, 1);
	E[6].set(V, E, F, 1, 6, 7, 5,14, 1);
	E[7].set(V, E, F, 6,10, 8, 6,39, 1);
	E[8].set(V, E, F,10, 5, 9, 7,30, 1);
	E[9].set(V, E, F, 5, 0, 5, 8,26, 1);
	
	E[10].set(V, E, F, 1, 2, 11, 14, 3,2);
	E[11].set(V, E, F, 2, 7, 12, 10,19,2);
	E[12].set(V, E, F, 7,11, 13, 11,44,2);
	E[13].set(V, E, F,11, 6, 14, 12,35,2);
	E[14].set(V, E, F, 6, 1, 10, 13, 6,2);
	
	E[15].set(V, E, F, 2, 3, 16, 19, 2,3);
	E[16].set(V, E, F, 3, 8, 17, 15,24,3);
	E[17].set(V, E, F, 8,12, 18, 16,49,3);
	E[18].set(V, E, F,12, 7, 19, 17,40,3);
	E[19].set(V, E, F, 7, 2, 15, 18,11,3);
	
	E[20].set(V, E, F, 3, 4, 21, 24, 1,4);
	E[21].set(V, E, F, 4, 9, 22, 20,29,4);
	E[22].set(V, E, F, 9,13, 23, 21,54,4);
	E[23].set(V, E, F,13, 8, 24, 22,45,4);
	E[24].set(V, E, F, 8, 3, 20, 23,16,4);
	
	E[25].set(V, E, F, 4, 0, 26, 29, 0,5);
	E[26].set(V, E, F, 0, 5, 27, 25, 9,5);
	E[27].set(V, E, F, 5,14, 28, 26,34,5);
	E[28].set(V, E, F,14, 9, 29, 27,50,5);
	E[29].set(V, E, F, 9, 4, 25, 28,21,5);


	
	E[30].set(V, E, F, 5,10, 31, 34, 8,6);
	E[31].set(V, E, F,10,15, 32, 30,38,6);
	E[32].set(V, E, F,15,19, 33, 31,59,6);
	E[33].set(V, E, F,19,14, 34, 32,51,6);
	E[34].set(V, E, F,14, 5, 30, 33,27,6);

	E[35].set(V, E, F, 6,11, 36, 39,13,7);
	E[36].set(V, E, F,11,16, 37, 35,43,7);
	E[37].set(V, E, F,16,15, 38, 36,55,7);
	E[38].set(V, E, F,15,10, 39, 37,31,7);
	E[39].set(V, E, F,10, 6, 35, 38, 7,7);

	E[40].set(V, E, F, 7,12, 41, 44,18,8);
	E[41].set(V, E, F,12,17, 42, 40,48,8);
	E[42].set(V, E, F,17,16, 43, 41,56,8);
	E[43].set(V, E, F,16,11, 44, 42,36,8);
	E[44].set(V, E, F,11, 7, 40, 43,12,8);

	E[45].set(V, E, F, 8,13, 46, 49,23,9);
	E[46].set(V, E, F,13,18, 47, 45,53,9);
	E[47].set(V, E, F,18,17, 48, 46,57,9);
	E[48].set(V, E, F,17,12, 49, 47,41,9);
	E[49].set(V, E, F,12, 8, 45, 48,17,9);

	E[50].set(V, E, F, 9,14, 51, 54,28,10);
	E[51].set(V, E, F,14,19, 52, 50,33,10);
	E[52].set(V, E, F,19,18, 53, 51,58,10);
	E[53].set(V, E, F,18,13, 54, 52,46,10);
	E[54].set(V, E, F,13, 9, 50, 53,22,10);



	E[55].set(V, E, F,15,16, 56, 59,37,11);
	E[56].set(V, E, F,16,17, 57, 55,42,11);
	E[57].set(V, E, F,17,18, 58, 56,47,11);
	E[58].set(V, E, F,18,19, 59, 57,52,11);
	E[59].set(V, E, F,19,15, 55, 58,32,11);

	F[0].from = &E[0];		F[0].Norm = Vec(0, 0, -1);
	F[1].from = &E[05];		F[1].Norm = Vec(-2.*V[3].X.x, -2.*V[3].X.y, -1.) / s5;
	F[2].from = &E[10]; 	F[2].Norm = Vec(-2.*V[4].X.x, -2.*V[4].X.y, -1.) / s5;
	F[3].from = &E[15];		F[3].Norm = Vec(-2.*V[0].X.x, -2.*V[0].X.y, -1.) / s5;
	F[4].from = &E[20];		F[4].Norm = Vec(-2.*V[1].X.x, -2.*V[1].X.y, -1.) / s5;
	F[5].from = &E[25];		F[5].Norm = Vec(-2.*V[2].X.x, -2.*V[2].X.y, -1.) / s5;
	F[6].from = &E[30];		F[6].Norm = -F[3].Norm;
	F[7].from = &E[35];		F[7].Norm = -F[4].Norm;
	F[8].from = &E[40];		F[8].Norm = -F[5].Norm;
	F[9].from = &E[45];		F[9].Norm = -F[1].Norm;
	F[10].from =&E[50];		F[10].Norm = -F[2].Norm;
	F[11].from =&E[55];		F[11].Norm = Vec(0, 0, 1.);

	for (int i=0; i<numV; i++)
		V[i].X *= .5;


	for (int i=0; i<numV; i++)
		V[i].X = Vec(V[i].X.x*Siz.x, V[i].X.y*Siz.y, V[i].X.z*Siz.z) * Ori + Pos;

	for (int f=0; f<numF; f++)
		F[f].update();

	for (int i=0; i<numV; i++)
		V[i].X *= 4/sqrt(10-sqrt(20));


	consistsOfOnlyTriangles = false;

	/*static int hej = 0;
	if (hej++ == 1)
		subdivide1();*/
}
