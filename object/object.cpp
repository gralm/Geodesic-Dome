#include "object.hpp"

using namespace std;

std::list<ObjectFN*> World::Objs;

void Edge::set(Vec* _V, Edge *_E, Face *_F, int _fr, int _to, int _next, int _prev, int _oppo, int _face)
{
	fr = _V + _fr;
	to = _V + _to;
	next = _E + _next;
	prev = _E + _prev;
	oppo = _E + _oppo;
	face = _F + _face;
}

ObjectFN::ObjectFN()
{
	V = 0;
	E = 0;
	F = 0;
	numV = numE = numF = 0;
}

ObjectFN::ObjectFN(int _vert, int _edge, int _face)
{
	numV = _vert;
	numE = _edge;
	numF = _face;
	V = new Vec[_vert];
	E = new Edge[_edge];
	F = new Face[_face];
}

ObjectFN::~ObjectFN()
{
	if (V && E && F)
	{
		cout << "ska döda obketet" << endl;
		delete[] V;
		delete[] E;
		delete[] F;
		cout << "Dödade objektet" << endl;
	} else {
		cout << "Går inte ta bort objekt, objektet är inte initierat propert " << endl;
	}
}

Face* ObjectFN::getFaces(int &numOfFaces) const
{
	numOfFaces = numF;
	return this->F;
}


void ObjectFN::print()
{
	cout << "antalV: " << numV << "\t";
	cout << "antalE: " << numE << "\t";
	cout << "antalF: " << numF << endl;
}


ObjCubeFN::ObjCubeFN(const Vec &Pos, const Vec &Siz, const Mat &Ori) :ObjectFN(8, 24, 6)
{
		cout << "Ska skapa en Kub" << endl;
	V[0] = Vec(-.5,		.5, 	-.5);
	V[1] = Vec(.5, 		.5,	 	-.5);
	V[2] = Vec(-.5,		.5, 	.5);
	V[3] = Vec(.5, 		.5, 	.5);
	V[4] = Vec(-.5,		-.5,	-.5);
	V[5] = Vec(.5, 		-.5,	-.5);
	V[6] = Vec(-.5,		-.5,	.5);
	V[7] = Vec(.5, 		-.5,	.5);


	for (int i=0; i<8; i++)
		V[i] = Vec(V[i].x*Siz.x, V[i].y*Siz.y, V[i].z*Siz.z) * Ori + Pos;

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
	cout << "Skapade en kub" << endl;

	for (int i=0; i<8; i++)
		cout << V[i] << endl;
	cout << "antalet killar \n";
	print();
}



ObjTetrahedronFN::ObjTetrahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori): ObjectFN(4, 12, 4)
{
	TYP s2 = static_cast<TYP>(sqrt(2.0));
	TYP s3 = static_cast<TYP>(sqrt(3.0));

	V[0] = Vec(1/s3, 		0, 		-0.5/(s3*s2));
	V[1] = Vec(-0.5/s3, 	0.5, 	-0.5/(s3*s2));
	V[2] = Vec(-0.5/s3, 	-0.5, 	-0.5/(s3*s2));
	V[3] = Vec(0, 			0, 		(s3*s2)*.25);

	for (int i=0; i<4; i++)
		V[i] = Vec(V[i].x*Siz.x, V[i].y*Siz.y, V[i].z*Siz.z) * Ori + Pos;


	E[0].set(V, E, F, 0, 1, 1, 2,11, 0);
	E[1].set(V, E, F, 1, 3, 2, 0, 8, 0);
	E[2].set(V, E, F, 3, 0, 0, 1, 5, 0);
	E[3].set(V, E, F, 3, 2, 4, 4, 7, 1);
	E[4].set(V, E, F, 2, 0, 5, 5, 9, 1);
	E[5].set(V, E, F, 0, 3, 3, 3, 2, 1);
	E[6].set(V, E, F, 1, 2, 7, 7,10, 2);
	E[7].set(V, E, F, 2, 3, 8, 8, 3, 2);
	E[8].set(V, E, F, 3, 1, 6, 6, 1, 2);
	E[9].set(V, E, F, 0, 2,10,11, 4, 3);
	E[10].set(V,E, F, 2, 1,11, 9, 6, 3);
	E[11].set(V,E, F, 1, 0, 9,10, 0, 3);

	F[0].from = &E[0];		F[0].Norm = (V[2] * (-2.0*s2/s3)) * Ori;
	F[1].from = &E[3];		F[1].Norm = (V[1] * (-2.0*s2/s3)) * Ori;
	F[2].from = &E[6];		F[2].Norm = (V[0] * (-2.0*s2/s3)) * Ori;
	F[3].from = &E[9];		F[3].Norm = (V[3] * (-2.0*s2/s3)) * Ori;
}


	//OBJ_DODECAHEDRON
ObjDodecahedronFN::ObjDodecahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori): ObjectFN(20, 60, 12)
{
	TYP s5 = static_cast<TYP>(sqrt(5.0));
	TYP sin1 = sqrt((5.0+s5) / 8.);
	TYP sin2 = sqrt((5.0-s5) / 8.);

		// A
	V[0] = Vec(1.0, 0.0, 0.0);
	V[1] = Vec(0.25*(s5-1.0), sin1, 0.0);
	V[2] = Vec(-0.25*(s5+1.0), sin2, 0.0);
	V[3] = Vec(V[2].x, -V[2].y, 0.0);
	V[4] = Vec(V[1].x, -V[1].y, 0.0);

		// B
	for (int i=0; i<5; i++)
		V[5+i] = V[i]*(0.5*(1+s5)) + Vec(0.0, 0.0, 0.25*(1.0-s5));
		//V[5+i] = V[i]*(0.25*(1+s5)) + Vec(0.0, 0.0, 0.25*(1.0-s5));

		// C
	for (int i=0; i<5; i++)
		V[10+i] = -V[i<2? i+8: i+3];

		// A
	for (int i=0; i<5; i++)
		V[i].z -= 0.25*(3+s5);

		// D
	for (int i=0; i<5; i++)
		V[15+i] = -V[(i<2)? i+3: i-2];


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


	
	E[30].set(V, E, F, 5,10, 31, 35, 8,6);
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
	F[1].from = &E[05];		F[1].Norm = Vec(-2.*V[3].x, -2.*V[3].y, -1.) / s5;
	F[2].from = &E[10]; 	F[2].Norm = Vec(-2.*V[4].x, -2.*V[4].y, -1.) / s5;
	F[3].from = &E[15];		F[3].Norm = Vec(-2.*V[0].x, -2.*V[0].y, -1.) / s5;
	F[4].from = &E[20];		F[4].Norm = Vec(-2.*V[1].x, -2.*V[1].y, -1.) / s5;
	F[5].from = &E[25];		F[5].Norm = Vec(-2.*V[2].x, -2.*V[2].y, -1.) / s5;
	F[6].from = &E[30];		F[6].Norm = -F[3].Norm;
	F[7].from = &E[35];		F[7].Norm = -F[4].Norm;
	F[8].from = &E[40];		F[8].Norm = -F[5].Norm;
	F[9].from = &E[45];		F[9].Norm = -F[1].Norm;
	F[10].from =&E[50];		F[10].Norm = -F[2].Norm;
	F[11].from =&E[55];		F[11].Norm = Vec(0, 0, 1.);

	for (int i=0; i<20; i++)
		V[i] *= .5;
}


//std::vector<ObjectFN*> Objs;
//public:
bool World::addObjectFN(int objType, const Vec &Pos, const Vec &Siz, const Mat &Ori)
{

	ObjectFN *nyFN = 0;
	switch(objType)
	{
		case OBJ_CUBE:{
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_TRUNCATED_ICOSAHEDRON:{
			cout << "OBJ_TRUNCATED_ICOSAHEDRON existerar inte ännu :(" << endl;
			break;
		}
		case OBJ_TETRAHEDRON: {
			nyFN = new ObjTetrahedronFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_DODECAHEDRON: {
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
			break;
		}
		default:
			return false;
	}

	if (nyFN == 0)
		return false;


	Objs.push_back(nyFN);

	cout << "Antal Objekt 1: " << Objs.size() << endl;

	return true;
}

bool World::removeAllObjects()
{
	int i=0;
	for (std::list<ObjectFN*>::iterator itO = Objs.begin(); itO != Objs.end(); itO++)
	{
		cout << "ska deleta obj num: " << i << endl;
		ObjectFN *elm = *itO;
		
		delete elm;
		cout << "deletat done" << endl;
	}

	return true;
}


	ObjectFN* World::getAnObject()
	{
		return static_cast<ObjectFN*>(*Objs.begin());
	}

