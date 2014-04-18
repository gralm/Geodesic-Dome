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
	E[10].set(V, E, F,2, 1,11, 9, 6, 3);
	E[11].set(V, E, F,1, 0, 9,10, 0, 3);

	F[0].from = &E[0];		F[0].Norm = (V[2] * (-2.0*s2/s3)) * Ori;
	F[1].from = &E[3];		F[1].Norm = (V[1] * (-2.0*s2/s3)) * Ori;
	F[2].from = &E[6];		F[2].Norm = (V[0] * (-2.0*s2/s3)) * Ori;
	F[3].from = &E[9];		F[3].Norm = (V[3] * (-2.0*s2/s3)) * Ori;
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
			break;
		}
		case OBJ_TETRAHEDRON: {
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
		cout << "Antal Objekt 2: " << Objs.size() << endl;
		return static_cast<ObjectFN*>(*Objs.begin());
	}

