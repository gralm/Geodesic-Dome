#include "object.hpp"
#include "objSimple.hpp"
#include "objCube.hpp"
#include "objTetrahedron.hpp"
#include "objDodecahedron.hpp"
#include "objTruncatedIcosahedron.hpp"
		//objTruncatedIcosahedron

#define ERR_PRINT(str_) 	cout << numOfErrs << ": " << str_ << endl; if (numOfErrs++ >= maxNumOfErrs) return false

using namespace std;

std::list<ObjectFN*> World::Objs;



int Vertex::countEdges()
{
	int n = 1;
	for (Edge *iterE = from->oppo->next; iterE != from && n<20; iterE = iterE->oppo->next)
	{
		n++; 
	}
	return n;
}


void Edge::set(Vertex* _V, Edge *_E, Face *_F, int _fr, int _to, int _next, int _prev, int _oppo, int _face)
{
	fr = _V + _fr;
	to = _V + _to;
	next = _E + _next;
	prev = _E + _prev;
	oppo = _E + _oppo;
	face = _F + _face;
	_F[_face].from = this;	//
	_V[_fr].from = this;
}

TYP Edge::length()
{
	Vec _l = to->X - fr->X;
	return sqrt(_l*_l);
}

void Face::update()
{
	Norm = (from->to->X - from->fr->X) & (from->next->to->X - from->next->fr->X);
	if (Norm * Norm > 0.00000001)
		Norm.norm();
	else 
		cout << "För dålig precision för att skapa normal i facet" << endl;
}


Vec Face::getCenter()
{
	int n = 0;
	Edge *iterE = from;
	Vec Ret_ = Vec(0,0,0);
	do {
		Ret_ += iterE->fr->X;
		n++;
		iterE = iterE->next;
	} while(iterE != from);

	return Ret_ / n;
}

int Face::countEdges()
{
	int k = 0;
	Edge *iterE = from;
	do {
		iterE = iterE->next;
		k++;
	} while(iterE != from);
	return k;
}

TYP Face::maxSinErr(int &N)		// N = num of edges
{
	TYP maxErr2 = 0.;
	N = 0;
	Edge *iterE = from;

	do {
		Vec thisVec = (iterE->to->X - iterE->fr->X);
		TYP val2 = thisVec * Norm;
		val2 = val2*val2 / (thisVec*thisVec);

		if (val2 > maxErr2)
			maxErr2 = val2;

		iterE = iterE->next;
		if (++N > 20)
			return -1.0;

	} while(iterE != from);

	return sqrt(maxErr2);
}

void Face::moveUp(TYP height)
{
	Edge *iterE = from;
	do {
		iterE->fr->X += Norm*height;
		iterE = iterE->next;
	} while(iterE != from);
}

void Face::rotate(TYP rad)
{
	double c = cos(rad);
	double s = sin(rad);
	Vec C_ = getCenter();
	Edge *iterE = from;
	do {
		Vec dNorm = Norm * (iterE->fr->X * Norm);
		iterE->fr->X = dNorm + (iterE->fr->X - dNorm)*c + (Norm & iterE->fr->X)*s;
		//iterE->fr->X = iterE->fr->X * c + (Norm & iterE->fr->X) * s;
		iterE = iterE->next;
	} while (iterE != from);
}


void Edge::print(Vertex *V0, Edge *E0, Face *F0)
{
	cout << "E[" << (this - E0) << "] {fr=" << (fr-V0) << ",\tto=" << (to-V0);
	cout << ",\tnext=" << (next-E0) << ",\tprev=" << (prev-E0) << ",\toppo=" << (oppo-E0);
	cout << ",\tface=" << (face-F0) << "}" << endl;
}


/*
Face *ObjectFN::truncatedFaceNum(Face *_Start, int _N, int _r, int _n)
{
	return 
}*/

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

	V = new Vertex[_vert];
	E = new Edge[_edge];
	F = new Face[_face];

	consistsOfOnlyTriangles = false;
}

ObjectFN::~ObjectFN()
{
	if (V && E && F)
	{
		//cout << "ska döda obketet" << endl;
		delete[] V;
		delete[] E;
		delete[] F;
		//cout << "Dödade objektet" << endl;
	} else {
		cout << "Går inte ta bort objekt, objektet är inte initierat propert " << endl;
	}
}

Face* ObjectFN::getFaces(int &numOfFaces) const
{
	numOfFaces = numF;
	return this->F;
}

Vec ObjectFN::getCenter()
{
	Vec Mitt = Vec(0, 0, 0);
	for (int v=0; v<numV; v++)
		Mitt += V[v].X;
	return Mitt / numV;
}

TYP ObjectFN::normalizeRadius()
{
	Vec Mitt = Vec(0, 0, 0);
	for (int v=0; v<numV; v++)
		Mitt += V[v].X;
	Mitt /= numV;

	for (int v=0; v<numV; v++)
		V[v].X -= Mitt;

	TYP rad2 = 0;
	for (int v=0; v<numV; v++)
		rad2 += V[v].X*V[v].X;
	TYP rad = sqrt(rad2/numV);

	cout << "Mitten: " << Mitt << endl;
	cout << "radie: " << rad << endl;

	for (int v=0; v<numV; v++) 
		V[v].X = V[v].X * (rad / sqrt(V[v].X * V[v].X)) + Mitt;

	for (int f=0; f<numF; f++)
		F[f].update();

	return rad;
}

void ObjectFN::edgeCompare(TYP &shortest, TYP &longest)
{
	Vec thisE = (E[0].to->X - E[0].fr->X);
	TYP longest2 = thisE*thisE;
	TYP shortest2 = thisE*thisE;
	for (int e=1; e<numE; e++)
	{
		Vec thisE = (E[e].to->X - E[e].fr->X);
		TYP thisE2 = thisE*thisE;
		longest2 = max(thisE2, longest2);
		shortest2 = min(thisE2, shortest2);
	}

	shortest = sqrt(shortest2);
	longest = sqrt(longest2);
	//return sqrt(longest2) - sqrt(shortest2);
}

bool ObjectFN::updateConsistsOfOnlyTriangles()
{

	consistsOfOnlyTriangles = true;

	for (int f=0; f<numF; f++)
	{
		int n = 0;
		Edge *itE = F[f].from;
		Edge *endE = itE->prev;

		while (itE != endE) {itE = itE->next; n++;}
		//cout << "n: " << n << "\tf:" << f << endl;
		if (n != 2)
			return consistsOfOnlyTriangles = false;
	}

	return consistsOfOnlyTriangles;
}

bool ObjectFN::CopyVEF(Vertex *nyV, Edge *nyE, Face *nyF)
{
	if (nyV) 
	{
		for (int v=0; v<numV; v++)
		{
			nyV[v].from = nyE + (V[v].from - E);
			nyV[v].X = V[v].X;
			nyV[v].Norm = V[v].Norm;
		}
	}

	if (nyE) 
	{
		for (int e=0; e<numE; e++)
		{
			nyE[e].fr = nyV + (E[e].fr - V);
			nyE[e].to = nyV + (E[e].to - V);
			nyE[e].next = nyE + (E[e].next - E);
			nyE[e].prev = nyE + (E[e].prev - E);
			nyE[e].oppo = nyE + (E[e].oppo - E);
			nyE[e].face = nyF + (E[e].face - F);
		}	
	}		

	if (nyF) 
	{
		for (int f=0; f<numF; f++)
		{
			nyF[f].from = nyE + (F[f].from - E);
			nyF[f].Norm = F[f].Norm;
		}
	}

	return true;
}

bool ObjectFN::EmptyVEF(Vertex *nyV, Edge *nyE, Face *nyF)
{
	for (int v=0; v<numV; v++)
		nyV[v].from = 0;

	for (int e=0; e<numE; e++)
	{
		nyE[e].fr = nyE[e].to = 0;
		nyE[e].next = nyE[e].prev = nyE[e].oppo = 0;
		nyE[e].face = 0;
	}

	for (int f=0; f<numF; f++)
		nyF[f].from = 0;

	return true;
}

bool ObjectFN::transform(const Vec *Pos, const Vec *Siz, const Mat *Ori)
{
	Vec preC_ = getCenter();
	Vec postC_ = (Pos? (preC_ + *Pos): preC_);
	
	for (int i=0; i<numV; i++)
	{
		V[i].X -= preC_;
		if (Siz)
			V[i].X = Vec(V[i].X.x * Siz->x, V[i].X.y * Siz->y, V[i].X.z * Siz->z);

		if (Ori)
			V[i].X = V[i].X * (*Ori);

		V[i].X += postC_;
	}
	return true;
} 

bool ObjectFN::normalizeNormals()
{
	Vec E1, E2;
	for (int f=0; f<numF; f++)
	{
		E1 = F[f].from->to->X - F[f].from->fr->X;
		E2 = F[f].from->next->to->X - F[f].from->next->fr->X;
		E1 = E1 & E2;
		E1.norm();
		F[f].Norm = E1;
	}
}

bool ObjectFN::test(unsigned int shapeType_) const
{

	int maxNumOfErrs = 30;
	int numOfErrs = 0;
	int maxNumEdgesPerFace = 10;
	TYP minEdgeLenSq = (E[0].to->X-E[0].fr->X) * (E[0].to->X-E[0].fr->X);
	TYP maxEdgeLenSq = minEdgeLenSq;
	int minEdgeLenId = 0;
	int maxEdgeLenId = 0;
	Vec _Center(0, 0, 0);

	int numFaceWithVerts[20];
	for (int i=0; i<20; i++)
		numFaceWithVerts[i] = 0;
	
	for (int v=0; v<numV; v++)
	{
		if (V[v].from < E || V[v].from >= E + numE)
			cout << "V[" << v << "].from outside legal interval" << endl;
		else if (V[v].from->fr != &V[v])
			cout << "V[" << v << "].from->fr != &V[" << v << "]" << endl;
		_Center += V[v].X;

	}
	_Center /= numV;



	for (int e=0; e<numE; e++)
	{
		TYP edgeLen = (E[e].to->X - E[e].fr->X) * (E[e].to->X - E[e].fr->X);

		/*if (edgeLen < minEdgeLenSq) {
			minEdgeLenSq = edgeLen;
			minEdgeLenId = e;
		}
		
		if (edgeLen < minEdgeLenSq) {
			minEdgeLenSq = edgeLen;
			maxEdgeLenId = e;
		}*/
		if (edgeLen < minEdgeLenSq) {
			minEdgeLenSq = edgeLen;
			minEdgeLenId = e;
		}
		
		if (edgeLen > maxEdgeLenSq) {
			maxEdgeLenSq = edgeLen;
			maxEdgeLenId = e;
		}

			// kolla om någon är utanför intervall.
		if (E[e].fr < V || E[e].fr >= V+numV){
			ERR_PRINT("E[" << e << "].fr outside legal interval");
		}
		
		if (E[e].to < V || E[e].to >= V+numV){
			ERR_PRINT("E[" << e << "].to outside legal interval");
		}

		if (E[e].next < E || E[e].next >= E+numE){
			ERR_PRINT("E[" << e << "].next outside legal interval");
		}

		if (E[e].prev < E || E[e].prev >= E+numE){
			ERR_PRINT("E[" << e << "].prev outside legal interval");
		}

		if (E[e].oppo < E || E[e].oppo >= E+numE){
			ERR_PRINT("E[" << e << "].oppo outside legal interval");
		}

		if (E[e].face < F || E[e].face >= F+numF){
			ERR_PRINT("E[" << e << "].face outside legal interval");
		}





		if (E[e].next->fr != E[e].to) {
			ERR_PRINT("E[" << e << "].next->fr != E[" << e << "].to");
		}

		if (E[e].prev->to != E[e].fr) {
			ERR_PRINT("E[" << e << "].prev->to != E[" << e << "].fr");
		}


		if (E[e].face != E[e].next->face) {
			ERR_PRINT("E[" << e << "].face != E[" << e << "].next->face");
		}

		if (E[e].next->prev != &E[e]) {
			ERR_PRINT("E[" << e << "].next->prev != &E[" << e << "]");
		}
		
		if (E[e].prev->next != &E[e]) {
			ERR_PRINT("E[" << e << "].prev->next != &E[" << e << "]");
		}
		
		if (E[e].oppo->oppo != &E[e]) {
			ERR_PRINT("E[" << e << "].oppo->oppo not correct");
		}

		for (int e2=0; e2<numE; e2++)
		{
			if (e2 == e)
				continue;

			if (E[e].fr == E[e2].fr && E[e].to == E[e2].to) {
				ERR_PRINT("E[" << e << "] connects same edge as E[" << e2 << "]");
			}

			if (E[e].fr == E[e2].to && E[e].to == E[e2].fr && E[e].oppo != &E[e2]) {
				ERR_PRINT("E["<< e << "] is connected as E[" << e2 << "] but is not each others opposite");
			}

			if (E[e].next == E[e2].next) {
				ERR_PRINT("E[" << e << "].next = E[" << e2 << "].next");
			}

			if (E[e].prev == E[e2].prev) {
				ERR_PRINT("E[" << e << "].prev = E[" << e2 << "].prev");
			}

			if (E[e].oppo == E[e2].oppo) {
				ERR_PRINT("E[" << e << "].oppo = E[" << e2 << "].oppo");
			}
		}

		
		Edge *iterStart = E[e].next;
		Edge *iter = &E[e];
		int i;
		for (i=0; i<maxNumEdgesPerFace && iter != iterStart; i++)
		{
			iter = iter->next;
		}

		if (i<1 || i>=maxNumEdgesPerFace-1) {
			ERR_PRINT("E[" << e << "].next1->next2->next3...->nextN = start, N = " << i);
		}
	}

	cout << "Longsta Edge: " << sqrt(maxEdgeLenSq) << "\t med len: " << maxEdgeLenId << endl;
	cout << "Minsta Edge: " << sqrt(minEdgeLenSq) << "\t med len: " << minEdgeLenId << endl;

		// check face normal
	for (int f=0; f<numF; f++)
	{
		Edge *itE = F[f].from;
		Edge *endE = itE->prev;
		Vec _FaceCenter = endE->fr->X;
		int i = 1;
		while (itE != endE) {
			if (itE->face != F + f)
			{
				ERR_PRINT("F[" << f << "].from har edge ansluten till F[" << (itE->face - F) << "]");
			}
			_FaceCenter += itE->fr->X;
			i++; 
			itE = itE->next;
		}

		numFaceWithVerts[i]++;

		if (shapeType_ & OBJ_TYPE_SPHERICAL) 
		{
			_FaceCenter /= i;
			_FaceCenter -= _Center;
			_FaceCenter.norm();

			TYP hej = _FaceCenter * F[f].Norm;
			if (hej < 0.7 || hej > 1.3)
			{
				cout << "Bad face in F[" << f << "]" << endl;
				cout << "\t F[" << f << "].Norm = " << F[f].Norm << endl;
				cout << "\t _FaceCenter = " << _FaceCenter << endl;
				cout << "\t hej = " << hej << endl << endl;
				if (numOfErrs++ >= maxNumOfErrs) return false;
			}
		}
	}

		// check that all faces has vertices on surface
	for (int f=0; f<numF; f++)
	{
		Edge *itE = F[f].from;
		Edge *endE = itE->prev;
		int i;
		while (itE != endE) {
			TYP _val = ((itE->to->X - itE->fr->X) * F[f].Norm) / sqrt((itE->to->X - itE->fr->X)*(itE->to->X - itE->fr->X));
			if (_val < -0.01 || _val > 0.01)
			{
				ERR_PRINT("F[" << f << "].Norm * E[" << (itE - E) << "] = " << _val);
			}
			i++; 
			itE = itE->next;
		}
	}

	for (int i=0; i<20; i++)
	{
		if (numFaceWithVerts[i])
			cout << "Antal Faces med " << i << " vertices = " << numFaceWithVerts[i] << endl;
	}
}

void ObjectFN::print()
{
	cout << "antalV: " << numV << "\t";
	cout << "antalE: " << numE << "\t";
	cout << "antalF: " << numF << endl;

	cout << "\tV[n] = {Edge *from, \tVec3 X, \tVec3 Norm}" << endl;
	for (int v=0; v<numV; v++)
		cout << "V[" << v << "] = {" << (V[v].from - E) << ",\t" << V[v].X << ",\t" << V[v].Norm << endl;

	cout << endl << "\tE[n] {Vertex *fr,\tVertex *to,\tEdge *next,\tEdge *prev,\tEdge *oppo,\tFace *face\tLength}" << endl;
	for (int e=0; e<numE; e++) {
		cout << "E[" << e << "] {" << (E[e].fr - V) << ",\t" << (E[e].to - V) << ",\t" << (E[e].next - E) << ",\t";
		cout << (E[e].prev - E) << ",\t" << (E[e].oppo - E) << ",\t" << (E[e].face - F) << ", " << E[e].length() << "}" << endl;
	}

	cout << endl << "\tF[n] {Edge *from,\tVec Norm}" << endl;
	for (int f=0; f<numF; f++)
	{
		Edge *iterE = F[f].from;
		int k=0;
		do {
			iterE = iterE->next;
			k++;
		} while(iterE != F[f].from);
		cout << "F[" << f << "] {" << F[f].from - E << ", " << F[f].Norm << ", " << k << "}" << endl;
	}

}


void ObjectFN::print(const Vertex *V_, const Edge *E_, const Face *F_, int numV_, int numE_, int numF_)
{
	cout << "antalV: " << numV_ << "\t";
	cout << "antalE: " << numE_ << "\t";
	cout << "antalF: " << numF_ << endl;

	cout << "\tV[n] = {Edge *from, \tVec3 X, \tVec3 Norm}" << endl;
	for (int v=0; v<numV_; v++)
		cout << "V[" << v << "] = {" << (V_[v].from - E_) << ",\t" << V_[v].X << ",\t" << V_[v].Norm << endl;

	cout << endl << "\tE[n] {Vertex *fr,\tVertex *to,\tEdge *next,\tEdge *prev,\tEdge *oppo,\tFace *face}" << endl;
	for (int e=0; e<numE_; e++) {
		cout << "E[" << e << "] {" << (E_[e].fr - V_) << ",\t" << (E_[e].to - V_) << ",\t" << (E_[e].next - E_) << ",\t";
		cout << (E_[e].prev - E_) << ",\t" << (E_[e].oppo - E_) << ",\t" << (E_[e].face - F_) << "}" << endl;
	}

	cout << endl << "\tF[n] {Edge *from,\tVec Norm, Edges}" << endl;
	for (int f=0; f<numF_; f++)
	{
		Edge *iterE = F_[f].from;
		int k=0;
		do {
			iterE = iterE->next;
			k++;
		} while(iterE != F_[f].from);
		cout << "F[" << f << "] {" << F_[f].from - E_ << ", " << F_[f].Norm << ", " << k << "}" << endl;
	}

}





bool ObjectFN::setPolygonHeight(TYP h, int N)
{
	Vec C_ = getCenter();
	cout << "Center: " << C_ << endl;

	for (int f=0; f<numF; f++)
	{
		if (F[f].countEdges() != N)
			continue;

		Edge *iterE = F[f].from;
		do {
			TYP presentHeight = (iterE->fr->X - C_)*F[f].Norm;
			iterE->fr->X += F[f].Norm * (h-presentHeight);
			iterE = iterE->next;
		} while(iterE != F[f].from);

	}
	return false;

}





	// bygger hörnen korrekt.
ObjectFN *ObjectFN::greenHousify(TYP b, TYP h)
{
	bool printar = false;
	Vec Ctr = getCenter();



	int numOf = numE/2;

	ObjectFN *Pin = new ObjectFN(12*numOf, 36*numOf, 8*numOf);
	numOf = 0;

	cout << "Tjena 01" << endl;
	for (int e=0; e<numE; e++)
	{
		if (E[e].fr > E[e].to)	
			continue;

		if (printar) 	cout << "E[" << e << "] = <" << (E[e].fr - V) << ", " << (E[e].to - V) << ">" << endl;
		Vec A1 = E[e].fr->X;

		Vec B1 = E[e].to->X;

			// A-vars
		Vec Z_A = (A1 - Ctr);

		Vec X_ = (B1 - A1);
		X_.norm();
		Vec Y_ = Z_A & X_;
		Y_.norm();

		Vec Z_ = X_ & Y_;
		Z_.norm();


			// (A11 - A1)*Y_ 	= b/2
			// (A11 - A1)*Y_1 	= -b/2
			// (A11 - A1)*Z_ 	= 0

			// 				[	Y_ 	]^T
			// (A11 - A1)	[	Y_1	]	=	b/2 [1 	-1 	0]
			//				[	Z_ 	]
			
			//								[	Y_1 x Z_	]
			// (A11 - A1) = b/2 [1 	-1 	0]	[	Z_ x Y_		] / (Y_ * (Y_1 x Z_))
			//								[	Y_ x Z_A	]

			// (A11 - A1) = b/2 (Z_ x (Y_ + Y_1)) / (Y_1*X_)
		Vec X_1 = (E[e].prev->fr->X - E[e].prev->to->X);
		//X_1.norm();
		Vec Y_1 = Z_A & X_1;
		Y_1.norm();
		Vec A11 = (Z_ & (Y_ + Y_1)) * (b / (2 * X_*Y_1));
		A11 += A1;


			// (A12 - A1)*Y_ 	= -b/2
			// (A12 - A1)*Y_2 	= b/2
			// (A12 - A1)*Z_ 	= 0
		Vec X_2 = (E[e].oppo->next->to->X - E[e].oppo->next->fr->X);
		//X_2.norm();
		Vec Y_2 = Z_A & X_2;
		Y_2.norm();
		Vec A12 = (Z_ & (Y_ + Y_2)) * (-b / (2 * X_*Y_2));
		A12 += A1;


		TYP alpha_ = h / (Z_A*Z_);

		Vec A2 = A1 - Z_A*alpha_;
		Vec A21 = A11 - Z_A*alpha_;
		Vec A22 = A12 - Z_A*alpha_;


			// B-vars
		Vec Z_B = (B1 - Ctr);

		X_1 = (E[e].next->to->X - E[e].next->fr->X);
		//X_1.norm();
		Y_1 = Z_B & X_1;
		Y_1.norm();
		Vec B11 = (Z_ & (Y_ - Y_1)) * (-b / (2 * X_*Y_1));
		B11 += B1;

		X_2 = (E[e].oppo->prev->fr->X - E[e].oppo->prev->to->X);
		//X_2.norm();
		Y_2 = Z_B & X_2;
		Y_2.norm();
		Vec B12 = (Z_ & (Y_2 - Y_)) * (-b / (2 * X_*Y_2));
		B12 += B1;

		alpha_ = h / (Z_B*Z_);

		Vec B2 = B1 - Z_B*alpha_;
		Vec B21 = B11-Z_B*alpha_;
		Vec B22 = B12-Z_B*alpha_;
		

			//////
		Pin->V[0 + 12*numOf].X = A1;		Pin->V[0 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[1 + 12*numOf].X = A12;		Pin->V[1 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[2 + 12*numOf].X = B12;		Pin->V[2 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[3 + 12*numOf].X = B1;		Pin->V[3 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[4 + 12*numOf].X = B11;		Pin->V[4 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[5 + 12*numOf].X = A11;		Pin->V[5 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[6 + 12*numOf].X = A2;		Pin->V[6 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[7 + 12*numOf].X = A22;		Pin->V[7 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[8 + 12*numOf].X = B22;		Pin->V[8 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[9 + 12*numOf].X = B2;		Pin->V[9 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[10 + 12*numOf].X = B21;		Pin->V[10 + 12*numOf].from = &Pin->E[0 + 36*numOf];
		Pin->V[11 + 12*numOf].X = A21;		Pin->V[11 + 12*numOf].from = &Pin->E[0 + 36*numOf];


		Pin->E[0 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 0, 1, 1, 5, 9, 0);
		Pin->E[1 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 1, 2, 2, 0, 13, 0);
		Pin->E[2 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 2, 3, 3, 1, 17, 0);
		Pin->E[3 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 3, 4, 4, 2, 21, 0);
		Pin->E[4 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 4, 5, 5, 3, 25, 0);
		Pin->E[5 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 5, 0, 0, 4, 29, 0);


		for (int i=0; i<6; i++)
		{
			int i4=i*4;
			int a0 = i;
			int a1 = 6+i;
			int a2 = 7+(i==5? -1: i);
			int a3 = 1+(i==5? -1: i);
			Pin->E[6+i4 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, a0, a1, 7+i4, 9+i4, (i==0? 28: i4+4),	1+i);
			Pin->E[7+i4 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, a1, a2, 8+i4, 6+i4, 35-i, 				1+i);
			Pin->E[8+i4 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, a2, a3, 9+i4, 7+i4, (i==5? 6: i4+10), 	1+i);
			Pin->E[9+i4 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, a3, a0, 6+i4, 8+i4, i, 				1+i);
		}


		Pin->E[30 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf,6, 11, 31, 35, 27, 7);
		Pin->E[31 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf,11,10, 32, 30, 23, 7);
		Pin->E[32 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf,10, 9, 33, 31, 19, 7);
		Pin->E[33 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 9, 8, 34, 32, 15, 7);
		Pin->E[34 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 8, 7, 35, 33, 11, 7);
		Pin->E[35 + 36*numOf].set(Pin->V + 12*numOf, Pin->E + 36*numOf, Pin->F + 8*numOf, 7, 6, 30, 34, 7, 7);

		numOf++;
	}

	for (int f=0; f<8*numOf; f++)
		Pin->F[f].update();

	return Pin;
}

bool ObjectFN::rotatePolygons(TYP angle, int N)			// Rotate faces with N vertices by angle
{
	int *inN = new int[numV];
	list<int> rotF;

	for (int v=0; v<numV; v++)
		inN[v] = 0;

	//if (printar)		cout << "numF = " << numF << endl;


		// kontrollera först att inga vertexes är bundna till två faces med N vertices
	for (int f=0; f<numF; f++)
	{
		int k=0;
		Edge *iterE = F[f].from;
		do {
			iterE = iterE->next;
			k++;
		} while(iterE != F[f].from);
		
		//cout << "F[" << f << "] = " << k << endl;
		if (k==N) {
			rotF.push_back(f);
			do {
				inN[iterE->fr - V]++;
				iterE = iterE->next;
			} while(iterE != F[f].from);
		}
	}

	for (int v=0; v<numV; v++) {
		//cout << "V[" << v << "] = " << inN[v] << endl;
		if (inN[v] > 1){
			delete[] inN;
			cout << "Cannot rotate polygons. V[" << v << "] connected to more than one face with " << N << "vertices" << endl;
			return false;
		}
	}

			// rotate

	TYP c = cos(angle);
	TYP s = sin(angle);
	//cout << "Center1: " << F[0].getCenter() << " = " << (V[21].X + V[22].X + V[24].X + V[26].X + V[28].X)*.2 << endl;


	cout << endl;


	for (list<int>::iterator iti = rotF.begin(); iti != rotF.end(); iti++)
	{
		Vec Cen = F[*iti].getCenter();

		Edge *iterE = F[*iti].from;
		do {
			iterE = iterE->next;
			iterE->fr->X = (iterE->fr->X - Cen)*c + ((iterE->fr->X - Cen) & F[*iti].Norm)*s + Cen;
		} while(iterE != F[*iti].from);
	}


	for (int f=0; f<numF; f++)
		F[f].update();

	delete[] inN;

	return true;
}

bool ObjectFN::snub(TYP h, TYP t, int N)	// rotera alla ytor med n kanter t radianer och expandera h
{
	int oldNumF = numF;
	expand(h);
	cout << "old num f = " << oldNumF << endl;
	cout << "ny num f = " << numF << endl;
	cout << "cos(t) = " << cos(t) << "\tsin(t) = " << sin(t) << endl;

	for (int n=0; n<oldNumF; n++)
	{
		if (F[n].countEdges() == N)
			F[n].rotate(t);
	}

	splitBrokenTetragons();
}

bool ObjectFN::splitBrokenTetragons()
{
	bool printar = false;
	bool *splitEdge = new bool[numF];

	int edges;
	int edgesToBeSplit = 0;
	for (int f=0; f<numF; f++)
	{
		splitEdge[f] = false;
		TYP err = F[f].maxSinErr(edges);
		if (err > 0.0001 && edges == 4) {
			splitEdge[f] = true;
			edgesToBeSplit++;
		}
		//cout << "F[" << f << "]:\t" << tjena << ",\t" << hej << ",\t" << F[f].countEdges() << endl;
	}

	Vertex *nyV = new Vertex[numV];
	Edge *nyE = new Edge[numE + 2*edgesToBeSplit];
	Face *nyF = new Face[numF + edgesToBeSplit];

	CopyVEF(nyV, nyE, nyF);
	int numFny = numF;

	for (int f=0; f<numF; f++)
	{
		if (splitEdge[f] == false)
			continue;


		Edge *iterE = nyF[f].from;

		Vec dist1 = iterE->next->to->X - iterE->fr->X;
		Vec dist2 = iterE->to->X - iterE->prev->fr->X;

		if (dist1*dist1 > dist2*dist2)
		{
			iterE = iterE->next;
			nyF[f].from = iterE;
		}


			// första sidan
		iterE->prev->next = &nyE[numE+1];
		iterE->prev->face = &nyF[numFny];
		iterE->prev->prev->prev = &nyE[numE+1];
		iterE->prev->prev->face = &nyF[numFny];
		nyE[numE+1].face = &nyF[numFny];
		nyF[numFny].from = &nyE[numE+1];
		nyE[numE+1].next = iterE->prev->prev;
		nyE[numE+1].prev = iterE->prev;
		nyE[numE+1].fr = iterE->fr;
		nyE[numE+1].to = iterE->next->to;
		nyE[numE+1].oppo = &nyE[numE];
		

			// andra sidan
		iterE->next->next = &nyE[numE];
		iterE->prev = &nyE[numE];
		nyE[numE].face = iterE->face;
		nyE[numE].next = iterE;
		nyE[numE].prev = iterE->next;
		nyE[numE].fr = iterE->next->to;
		nyE[numE].to = iterE->fr;
		nyE[numE].oppo = &nyE[numE+1];


		nyF[numFny].update();
		nyF[f].update();
		numFny++;
		numE += 2;
	}

	delete[] V;
	delete[] E;
	delete[] F;

	V = nyV;
	E = nyE;
	F = nyF;
	numF = numFny;
}



bool ObjectFN::alternate(bool removeFirstVertice)
{
		// kolla så att polyhedron är even sided
	int newNumOfFaces = 0;
	for (int f=0; f<numF; f++)
	{
		int numOfFaceEdges = F[f].countEdges();
		if (numOfFaceEdges < 4) {
			cout << "Inte possible att get alternation av polyhedron med mindre than 4 faces" << endl;
			return false;
		} else if (numOfFaceEdges&1) {
			cout << "Inte possible att get alternation av polyhedron med udda antal edges" << endl;
			return false;
		} else if (numOfFaceEdges > 4) {
			newNumOfFaces++;
		}
	}

	if (newNumOfFaces == 0)
	{
		cout << "Inte possible att get alternation av polyhedron utan edges med mer than 4 " << endl;
		return false;
	}
}


ObjectFN* World::addObjectFN(int objType, const Vec &Pos, const Vec &Siz, const Mat &Ori)
{

	ObjectFN *nyFN = 0;
	switch(objType)
	{
		case OBJ_SIMPLE:{
			nyFN = new ObjSimpleFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_CUBE:{
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_TRUNCATED_ICOSAHEDRON: {
			nyFN = new ObjTruncatedIcosahedronFN(Pos, Siz*sqrt(.5 + .5/sqrt(5.)), Ori);
			break;
		}
		case OBJ_TETRAHEDRON: {
			nyFN = new ObjTetrahedronFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_ICOSAHEDRON: {
				// fixa Siz-multiplikatorn här
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
			nyFN->makeDual();
			break;
		}
		case OBJ_DODECAHEDRON: {
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_TRUNCATED_TETRAHEDRON: {
			nyFN = new ObjTetrahedronFN(Pos, Siz*2, Ori);
			nyFN->subdivide1();
			nyFN->makeDual();
			break;	
		}
		case OBJ_TRUNCATED_OCTAHEDRON: {
			nyFN = new ObjCubeFN(Pos, Siz * sqrt(8./3), Ori);
			nyFN->subdivide1();
    		nyFN->makeDual();
    		break;
		}
		case OBJ_OCTAHEDRON: {
			nyFN = new ObjCubeFN(Pos, Siz * sqrt(9./8.), Ori);
			nyFN->makeDual();
			break;
		}
		case OBJ_TRUNCATED_CUBE: {
			nyFN = new ObjCubeFN(Pos, Siz * (1+sqrt(2)), Ori);
			nyFN->truncate(2. - sqrt(2));
			cout << "trunkerad" << endl;
			break;
		}
		case OBJ_ICOSIDODECAHEDRON: {
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
			
			nyFN->rectify();
			break;
		}
		case OBJ_CUBOCTAHEDRON: {
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			nyFN->rectify();
			break;
		}
		case OBJ_RHOMBICOSIDODECAHEDRON: {
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
			nyFN->expand(1/sqrt(2-2/sqrt(5)));
			break;
		}
		case OBJ_SNUBDODECAHEDRON_CW: 
		case OBJ_SNUBDODECAHEDRON_CCW: {
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
    		nyFN->expand(Siz.x*1.0);
			nyFN->rotatePolygons(.22874989202202764 * (objType==OBJ_SNUBDODECAHEDRON_CW? 1.0: -1.0), 5);
    		/*
    		if (objType == OBJ_SNUBDODECAHEDRON_CCW)
    			nyFN->rotatePolygons(-.22874989202202764, 5);
    		else
    			nyFN->rotatePolygons(.22874989202202764, 5);*/
    		nyFN->splitBrokenTetragons();
    		nyFN->setPolygonHeight(1.9809159472818407 * Siz.x, 5);
    		nyFN->normalizeNormals();
			break;
		}

		case OBJ_RHOMBICUBOCTAHEDRON: {
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			nyFN->expand(Siz.x*sqrt(.5));
			break;
		}

		case OBJ_GOLDBERG_1_2:
		{
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
    		nyFN->expand(Siz.x*1.0);
			nyFN->rotatePolygons(.22874989202202764, 5);
    		nyFN->splitBrokenTetragons();
    		nyFN->setPolygonHeight(1.9809159472818407 * Siz.x, 5);
    		nyFN->normalizeNormals();
    		//nyFN->subdivide1(5, 0.056440253827226845*2*Siz.x);
    		nyFN->subdivide1(5, 0.11288050765445369*Siz.x);
    		nyFN->makeDual();
			break;
		}

		case OBJ_SNUBCUBE_CW:
		case OBJ_SNUBCUBE_CCW:
		{
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			nyFN->snub(0.64261350892596258*Siz.x, (objType==OBJ_SNUBCUBE_CCW? 1.0: -1.0)*0.287413148757777626, 4);
			break;
		}

			// hexagons har vinklarna acos(1/sqr(3) ~= 54,7 deg och
		case OBJ_CHAMFERED_CUBE:
		{
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			nyFN->chamfer(Siz.x * 2/sqrt(3));
			break;
		}

		case OBJ_CHAMFERED_DODECAHEDRON:
		{
			//nyFN = World::addObjectFN(OBJ_DODECAHEDRON, Vec(0, 0, 0), Vec(1, 1, 1)*0.5, Mat(1,0,0, 0,1,0, 0,0,1));
			nyFN = new ObjDodecahedronFN(Pos, Siz, Ori);
    		//nyFN->chamfer(1.6180339888*Siz.x);
    		nyFN->chamfer(1.6180339887498927*Siz.x);
    		
    		break;
    	}


		default:
			cout << "finns inget like this objekt att plocka forward" << endl;
			return 0;
	}

	if (nyFN != 0)
		Objs.push_back(nyFN);

	return nyFN;
}

bool World::addObjectFN(ObjectFN *pObj)
{
	Objs.push_back(pObj);
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

std::list<ObjectFN*>* World::getObjectListPointer()
{
	return &Objs;
}

void ObjectFN::tabort()
{
	cout << "Face1: " << F[0].Norm << endl;
	cout << "Face2: " << F[0].from->oppo->face->Norm << endl;
	cout << "Face1*Face2: " << F[0].Norm * F[0].from->oppo->face->Norm << endl;
}