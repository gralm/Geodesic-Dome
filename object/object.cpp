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


Vertex *frV__ = 0;
Vertex *toV__ = 0;
Edge *frE__ = 0;
Edge *toE__ = 0;
Face *frF__ = 0;
Face *toF__ = 0;

void check(Vertex *v, int l)
{
	if (v < frV__ || v >= toV__)
		cout << l << ", Kass vertex: V[" << (v-frV__) << "]" << endl;
}

void check(Edge *e, int l)
{
	if (e < frE__ || e >= toE__)
		cout << l << ", Kass edge: E[" << (e-frE__) << "]" << endl;
}

void check(Face *f, int l)
{
	if (f < frF__ || f >= toF__)
		cout << l << ", Kass face: F[" << (f-frF__) << "]" << endl;
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
	int n = 1;
	Edge *iterE = from;
	Vec Ret_ = iterE->fr->X;
	do {
		iterE = iterE->next;
		Ret_ += iterE->fr->X;
		n++;
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


void Edge::print(Vertex *V0, Edge *E0, Face *F0)
{
	cout << "E[" << (this - E0) << "] {fr=" << (fr-V0) << ",\tto=" << (to-V0);
	cout << ",\tnext=" << (next-E0) << ",\tprev=" << (prev-E0) << ",\toppo=" << (oppo-E0);
	cout << ",\tface=" << (face-F0) << "}" << endl;
}


int ObjectFN::truncatedEdgeNum(int _N, int _r, int _n, int _p)
{
	int _ERR = -10000;
	int returnval = _ERR;

	if (_p<0 || _p>=3)
		returnval = _ERR;
	else if (_n<0 || _n>2*_r)
		returnval = _ERR;
	else if (_r >= _N)
		returnval = _ERR;
	else if (_r == _N-1) {
		if (_n&1) {
			returnval = 3*_r*_r - 2*_r + (5*_n - 3)/2 + _p;
		} else {
			if ((_p == 1) || (_n==0 && _p==0) || (_n==2*_r && _p==2))
				returnval = _ERR;
			else
				returnval = 3*_r*_r - 2*_r + (5*_n - 2 + _p)/2;
		}
		
	} else if (_r == 0)
		returnval = (_n || _p != 1)? (0): _ERR;
	else if (_r < 0)
		returnval = _ERR;
	else 
		returnval = 3*_r*_r - 2*_r + 3*_n - 1 + _p; 

		// felhantering
	//if (returnval == _ERR)
	//	cout << "Här räknades det fel Edge när _r = " << _r << ", _n = " << _n << ", _p = " << _p << endl;

	return returnval;
}


Edge *ObjectFN::truncatedEdgeNum(Edge *_Start, int _N, int _r, int _n, int _p)
{
	int _ERR = -1;
	int returnval = _ERR;

	if (_p<0 || _p>2 || _r<0 || _r>_N-1 || _n<0 || _n>2*_N)
		returnval = _ERR;
	else if ((_n == 0 && _p == 0) || (_n == 2*_r && _p == 2) || (_r == _N-1 && !(_n&1) && _p==1))
		returnval = _ERR;
	else if (_r < _N-1)
		returnval = 3*_r*_r - 2*_r + 3*_n - 1 + _p;
	else if (_n&1)
		returnval = 3*_r*_r - 2*_r + (5*_n-3)/2 + _p;
	else if (_p == 1)
		returnval = _ERR;
	else
		returnval = 3*_r*_r - 2*_r + (5*_n-2)/2 + _p/2;

	//if (returnval != _ERR)
	//cout << "N=" << _N << ", _r=" << _r << ", _n=" << _n << ", _p=" << _p << ", val = " << returnval << endl;

	return (returnval == _ERR)? 0: _Start+returnval;

}

Vertex *ObjectFN::truncatedVertexNum(Vertex *_Start, int _N, int _r, int _n, int _p)
{	
	int _ERR = -1;
	int returnval = _ERR;

	if (_n<0 || _n>_r*2)
		returnval = _ERR;
	else if (_r<1 || _r>=_N)
		returnval = _ERR;
	else if ((_r == _N-1) && ((_n&1 && _p==1) || (!(_n&1) && (_p != 2))))
		returnval = _ERR;
	else if (_n&1) {
		switch(_p) {
			case 0:
				returnval = (_r<1 || (_n<3))? _ERR: (_r-2)*(_r-1)/2 + (_n-3)/2;
				break;
			case 1:
				returnval = (_r<=0 && _r>=_N-1)? _ERR: (_r-1)*_r/2 + (_n-1)/2;
				break;
			case 2:
				returnval = (_n >= 2*_r-1)? _ERR: (_r-2)*(_r-1)/2 + (_n-1)/2;
				break;
			default:
				returnval = _ERR;
		}
	} else {
		switch(_p) {
			case 0:
				returnval = (_n==0)? _ERR: (_r-1)*_r/2 + _n/2 - 1;
				break;
			case 1:
				returnval = (_n == 2*_r)? _ERR: (_r-1)*_r/2 + _n/2;
				break;
			case 2:
				returnval = (_n == 0 || _n == 2*_r)? _ERR: (_r-2)*(_r-1)/2 + _n/2 - 1;
				break;
			default:
				returnval = _ERR;
		}

	}

	return (returnval == _ERR)? 0: _Start+returnval;
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

	cout << endl << "\tE[n] {Vertex *fr,\tVertex *to,\tEdge *next,\tEdge *prev,\tEdge *oppo,\tFace *face}" << endl;
	for (int e=0; e<numE; e++) {
		cout << "E[" << e << "] {" << (E[e].fr - V) << ",\t" << (E[e].to - V) << ",\t" << (E[e].next - E) << ",\t";
		cout << (E[e].prev - E) << ",\t" << (E[e].oppo - E) << ",\t" << (E[e].face - F) << "}" << endl;
	}

	cout << endl << "\tF[n] {Edge *from,\tVec Norm}" << endl;
	for (int f=0; f<numF; f++)
		cout << "F[" << f << "] {" << F[f].from - E << ", " << F[f].Norm << "}" << endl;

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

	cout << endl << "\tF[n] {Edge *from,\tVec Norm}" << endl;
	for (int f=0; f<numF_; f++)
		cout << "F[" << f << "] {" << F_[f].from - E_ << ", " << F_[f].Norm << "}" << endl;

}





	// skapa en vertex i mitten på varje yta. 
	// Justera den nya vertexen på rätt avstånd från mittpunkten
	// Skapa nya ytor faces 
bool ObjectFN::subdivide1()
{
	Vertex *nyV = new Vertex[numV + numF];
	Edge *nyE = new Edge[numE*3];
	Face *nyF = new Face[numE];

	//cout << "allokerar: " << (numV + numF) << " V \t";
	//cout << (numE*3) << " E \t";
	//cout << numE << " F" << endl;

	for (int i=0; i<numV; i++) {
		nyV[i].X = V[i].X;
		nyV[i].Norm = V[i].Norm;
		nyV[i].from =  nyE + (V[i].from - E);
	}
	
	for (int i=0; i<numE; i++) {
		nyE[i].fr = nyV + (E[i].fr - V);
		nyE[i].to = nyV + (E[i].to - V);
		nyE[i].next = nyE + (E[i].next - E);
		nyE[i].prev = nyE + (E[i].prev - E);
		nyE[i].oppo = nyE + (E[i].oppo - E);
		nyE[i].face = nyF + (E[i].face - F);
	}

	int numVny = numV;
	int numEny = numE;
	int numFny = 0;

	Vec M_ = Vec(0, 0, 0);
	for (int v=0; v<numV; v++)
		M_ += V[v].X;

	M_ /= static_cast<TYP>(numV);
	//cout << "Mitten: " << M_ << endl;

	for (int f=0; f<numF; f++)
	{
		Vec VfaceCenter_(0, 0, 0);
		////cout << "\tface: " << f << endl;
		Edge *iterE = F[f].from;
		int numFaceVerts = 0;
		TYP varians = 0;
		do {
			varians += (iterE->fr->X - M_)*(iterE->fr->X - M_);
			VfaceCenter_ += iterE->fr->X;
			//cout << (iterE - E) << "\t";
			iterE = iterE->next;
			numFaceVerts++;
		} while (iterE != F[f].from);
		//cout << endl;

		//cout << "numV: " << numVny << endl;
		//cout << "numE: " << numEny << endl;
		//cout << "numF: " << numFny << endl;

		varians /= static_cast<TYP>(numFaceVerts);
		//cout << " numFaceVerts: " << numFaceVerts << endl;
		//cout << " vad blir fel here? " << varians << endl;
		//cout <<  "snitt avstånd från centrum: " << sqrt(varians) << endl;

		VfaceCenter_ /= numFaceVerts;
		VfaceCenter_ -= M_;
		VfaceCenter_.norm();
		VfaceCenter_ *= sqrt(varians);
		VfaceCenter_ += M_;
		//cout << "Face center: " << VfaceCenter_ << endl;

			
		nyV[numVny].X = VfaceCenter_;
		nyV[numVny].from = 0;


		
			/// Fixa de nya fäjsen här:
		Edge *finalE = iterE = nyE + (F[f].from - E);
		Edge *prevFinalE = finalE->prev;

		do {
				// Skapa ny face i varje iteration
			nyF[numFny].from = iterE;

				// kom ihåg vilken som är nästa 
			Edge *nextE = iterE->next;
			iterE->next = &nyE[numEny];
			iterE->prev = &nyE[numEny+1];
			iterE->face = &nyF[numFny];

			//cout << "iterE[" << (iterE - nyE) << "]: " << "{ " << (iterE->fr - nyV) << ", " << (iterE->to - nyV) << ", " << (iterE->next - nyE)
			//		 << ", " << (iterE->prev - nyE) << ", " << (iterE->oppo - nyE) << ", " << (iterE->face - nyF) << "}" << endl;


			nyE[numEny].fr = iterE->to;
			nyE[numEny].to = &nyV[numVny];
			nyE[numEny].next = &nyE[numEny+1];
			nyE[numEny].prev = iterE;
			nyE[numEny].oppo = &nyE[numEny+3];
			nyE[numEny].face = &nyF[numFny];

			//cout << "iterE[" << (&nyE[numEny] - nyE) << "]: " << "{ " << (nyE[numEny].fr - nyV) << ", " << (nyE[numEny].to - nyV) << ", " << (nyE[numEny].next - nyE)
			//		 << ", " << (nyE[numEny].prev - nyE) << ", " << (nyE[numEny].oppo - nyE) << ", " << (nyE[numEny].face - nyF) << "}" << endl;


			nyE[numEny+1].fr = &nyV[numVny];
			nyE[numEny+1].to = iterE->fr;
			nyE[numEny+1].next = iterE;
			nyE[numEny+1].prev = &nyE[numEny];
			nyE[numEny+1].oppo = &nyE[numEny-2];
			nyE[numEny+1].face = &nyF[numFny];

			//cout << "iterE[" << (&nyE[numEny + 1] - nyE) << "]: " << "{ " << (nyE[numEny + 1].fr - nyV) << ", " << (nyE[numEny + 1].to - nyV) << ", " << (nyE[numEny + 1].next - nyE)
			//		 << ", " << (nyE[numEny + 1].prev - nyE) << ", " << (nyE[numEny + 1].oppo - nyE) << ", " << (nyE[numEny + 1].face - nyF) << "}" << endl;


			iterE = nextE;
			numEny += 2;
			numFny += 1;
		} while (iterE != finalE);

		iterE->prev->oppo = prevFinalE->next;
		prevFinalE->next->oppo = iterE->prev;

		//cout << "E[" << (iterE->prev - nyE) << "].oppo: " << (iterE->prev->oppo - nyE) << endl;
		//cout << "E[" << (prevFinalE->next - nyE) << "].oppo: " << (prevFinalE->next->oppo - nyE) << endl;

		//cout << "Färdig med denna facet" << endl;
		numVny++;
	}

	for (int e=0; e<numEny; e++) {
		nyE[e].fr->from = &nyE[e];
	}

	/*for (int v=0; v<numVny; v++) {
		cout << "V[" << v << "].from = E[" << nyV[v].from - nyE << "]" << endl;
	}*/


	delete[] V;
	delete[] E;
	delete[] F;

	V = nyV;
	E = nyE;
	F = nyF;

	numV = numVny;
	numE = numEny;
	numF = numFny;


	for (int f=0; f<numF; f++)
	{
		F[f].Norm = (F[f].from->to->X - F[f].from->fr->X) & (F[f].from->next->to->X - F[f].from->next->fr->X);
		F[f].Norm.norm();
		//cout << "Face[" << f << "]: " << F[f].Norm << endl;
	}
	consistsOfOnlyTriangles = true;
	return true;
}







	// http://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface. 
	// Catmull Clark Subdivision av ytor
bool ObjectFN::subdivide2()
{
	if (consistsOfOnlyTriangles == false)
	{
		cout << "Kan inte utföra subdivide2 innan alla faces är trianglar. Kör subdivide1 först" << endl;
		return false;
	}

	bool printar = false;

		// skapa nya vertex, edges och faces som ska användas.
	int numVny = numV;
	int numEny = numE;
	int numFny = numF;

	Vertex *nyV = new Vertex[numV + numE/2];
	Edge *nyE = new Edge[4*numE];
	Face *nyF = new Face[4*numF];

	for (int i=0; i<numV; i++) {
		nyV[i].X = V[i].X;
		nyV[i].Norm = V[i].Norm;
		nyV[i].from =  nyE + (V[i].from - E);
	}
	
	for (int i=0; i<numE; i++) {
		nyE[i].fr = nyV + (E[i].fr - V);
		nyE[i].to = nyV + (E[i].to - V);
		nyE[i].next = nyE + (E[i].next - E);
		nyE[i].prev = nyE + (E[i].prev - E);
		nyE[i].oppo = nyE + (E[i].oppo - E);
		nyE[i].face = nyF + (E[i].face - F);
	}

	for (int i=0; i<numF; i++)
		nyF[i].from = nyE + (F[i].from - E);

	if (printar) {
		cout << "*** innan vi sätter igång så ser det ut på följande vis ***" << endl;
			////////////////////////////////////////////////
		cout << "antalV: " << numV << "\t";
		cout << "antalE: " << numE << "\t";
		cout << "antalF: " << numF << endl;


		cout << "\tV[n] = {Edge *from, \tVec3 X, \tVec3 Norm}" << endl;
		for (int v_=0; v_<numVny; v_++)
			cout << "nyV[" << v_ << "] = {" << (nyV[v_].from - nyE) << ",\t" << nyV[v_].X << ",\t" << nyV[v_].Norm << endl;

		cout << endl << "\tE[n] {Vertex *fr,\tVertex *to,\tEdge *next,\tEdge *prev,\tEdge *oppo,\tFace *face}" << endl;
		for (int e_=0; e_<numEny; e_++) {
			cout << "nyE[" << e_ << "] {" << (nyE[e_].fr - nyV) << ",\t" << (nyE[e_].to - nyV) << ",\t" << (nyE[e_].next - nyE) << ",\t";
			cout << (nyE[e_].prev - nyE) << ",\t" << (nyE[e_].oppo?(nyE[e_].oppo - nyE) :-1) << ",\t" << (nyE[e_].face - nyF) << "}" << endl;
		}

		cout << endl << "\tF[n] {Edge *from,\tVec Norm}" << endl;
		for (int f_=0; f_<numFny; f_++)
			cout << "nyF[" << f_ << "] {" << nyF[f_].from - nyE << ", " << nyF[f_].Norm << "}" << endl;
	}
		/////////////////////////////////////////////////



	/////////////////////////////////////////////////////
	for (int f=0; f<numF; f++)
	//for (int f=0; f<3; f++)
	{


		if (printar) cout << endl << endl << "*********** Nytt face: " << f << endl;
		Edge *iterE = nyF[f].from;

		
		bool dividedOppo[3];
		Vertex *nyV_[3];

		for (int i=0; i<3; i++) {
			if ((iterE->fr==iterE->oppo->to) && (iterE->to==iterE->oppo->fr))
			{
				nyV_[i] = &nyV[numVny];
				nyV_[i]->X = (iterE->fr->X + iterE->to->X) * .5;
				iterE->to = nyV_[i];
				dividedOppo[i] = false;
				numVny++;
			} else {
				dividedOppo[i] = true;
				nyV_[i] = iterE->oppo->to;
			}
			if (printar) cout << "nyV_[" << i << "] = " << nyV_[i] - nyV << endl;

			if (printar) cout << "Edge[" << i << "] is " << (dividedOppo[i]? "divided": "undivided") << endl;
			iterE = iterE->next;
			//cout << "if ((" << (iterE->fr )
		}

		if (printar) cout << "numEny = " << numEny << " edges" << endl;
		if (printar) cout << "iterE = nyE[" << iterE - nyE << "]" << endl;

			// befintligt face
		nyE[numEny + 0].fr = nyV_[0];
		nyE[numEny + 0].to = nyV_[2];
		nyE[numEny + 0].next = &nyE[numEny + 1];
		nyE[numEny + 0].prev = iterE;
		nyE[numEny + 0].oppo = &nyE[numEny + 6];
		nyE[numEny + 0].face = &nyF[f];

		nyE[numEny + 1].fr = nyV_[2];
		nyE[numEny + 1].to = iterE->fr;	// iterE-> fr = 
		
		nyE[numEny + 1].next = iterE;
		nyE[numEny + 1].prev = &nyE[numEny + 0];
		if (dividedOppo[0]) {
			nyE[numEny + 3].oppo = iterE->oppo;
			iterE->oppo = nyE[numEny + 3].oppo->next->oppo->next->oppo->next;
			nyE[numEny + 3].oppo->oppo = &nyE[numEny + 3];
			iterE->oppo->oppo = iterE;
			if (printar) cout << "dividedoppo[0]" << endl;
			if (printar) cout << "nyE[" << numEny + 3 << "].oppo = nyE[" << nyE[numEny + 3].oppo - nyE << "]" << endl;
			if (printar) cout << "nyE[" << iterE - nyE << "].oppo = " << iterE->oppo - nyE << endl;
		} else {
			nyE[numEny + 5].oppo = 0;
		}
		nyE[numEny + 1].face = &nyF[f];

			// första nya facet
		nyE[numEny + 2].fr = nyV_[1];
		nyE[numEny + 2].to = nyV_[0];
		nyE[numEny + 2].next = &nyE[numEny + 3];
		nyE[numEny + 2].prev = iterE->next;
		nyE[numEny + 2].oppo = &nyE[numEny + 7];
		nyE[numEny + 2].face = &nyF[numFny + 0];

		nyE[numEny + 3].fr = nyV_[0];
		nyE[numEny + 3].to = iterE->next->fr;
		nyE[numEny + 3].next = iterE->next;
		nyE[numEny + 3].prev = &nyE[numEny + 2];
		if (dividedOppo[1]) {
			nyE[numEny + 5].oppo = iterE->next->oppo;
			iterE->next->oppo = nyE[numEny + 5].oppo->next->oppo->next->oppo->next;
			nyE[numEny + 5].oppo->oppo = &nyE[numEny + 5];
			iterE->next->oppo->oppo = iterE->next;
			if (printar) cout << "dividedoppo[1]" << endl;
			if (printar) cout << "nyE[" << numEny + 5 << "].oppo = nyE[" << nyE[numEny + 5].oppo - nyE << "]" << endl;
			if (printar) cout << "nyE[" << iterE->next - nyE << "].oppo = " << iterE->next->oppo - nyE << endl;
		} else {
			nyE[numEny + 5].oppo = 0;
		}
		nyE[numEny + 3].face = &nyF[numFny + 0];

		nyF[numFny].from = &nyE[numEny + 2];


			// andra nya facet
		nyE[numEny + 4].fr = nyV_[2];
		nyE[numEny + 4].to = nyV_[1];
		nyE[numEny + 4].next = &nyE[numEny + 5];
		nyE[numEny + 4].prev = iterE->prev;
		nyE[numEny + 4].oppo = &nyE[numEny + 8];
		nyE[numEny + 4].face = &nyF[numFny + 1];

		nyE[numEny + 5].fr = nyV_[1];
		nyE[numEny + 5].to = iterE->prev->fr;
		if (printar) cout << "HALLÅÅ: nyE[" << numEny + 5 << "] = " << nyE[numEny + 5].to - nyV << endl;
		nyE[numEny + 5].next = iterE->prev;
		nyE[numEny + 5].prev = &nyE[numEny + 4];
		if (dividedOppo[2]) {
			nyE[numEny + 1].oppo = iterE->prev->oppo;
			iterE->prev->oppo = nyE[numEny + 1].oppo->next->oppo->next->oppo->next;
			nyE[numEny + 1].oppo->oppo = &nyE[numEny + 1];
			iterE->prev->oppo->oppo = iterE->prev;
			if (printar) cout << "dividedoppo[2]" << endl;
			if (printar) cout << "nyE[" << numEny + 1 << "].oppo = nyE[" << nyE[numEny + 1].oppo - nyE << "]" << endl;
			if (printar) cout << "nyE[" << iterE->prev - nyE << "].oppo = " << iterE->prev->oppo - nyE << endl;
		} else {
			nyE[numEny + 5].oppo = 0;
		}
		nyE[numEny + 5].face = &nyF[numFny + 1];

		nyF[numFny + 1].from = &nyE[numEny + 4];



			// tredje nya facet
		nyE[numEny + 6].fr = nyV_[2];
		nyE[numEny + 6].to = nyV_[0];
		nyE[numEny + 6].next = &nyE[numEny + 7];
		nyE[numEny + 6].prev = &nyE[numEny + 8];
		nyE[numEny + 6].oppo = &nyE[numEny + 0];
		nyE[numEny + 6].face = &nyF[numFny + 2];

		nyE[numEny + 7].fr = nyV_[0];
		nyE[numEny + 7].to = nyV_[1];
		nyE[numEny + 7].next = &nyE[numEny + 8];
		nyE[numEny + 7].prev = &nyE[numEny + 6];
		nyE[numEny + 7].oppo = &nyE[numEny + 2];
		nyE[numEny + 7].face = &nyF[numFny + 2];

		nyE[numEny + 8].fr = nyV_[1];
		nyE[numEny + 8].to = nyV_[2];
		nyE[numEny + 8].next = &nyE[numEny + 6];
		nyE[numEny + 8].prev = &nyE[numEny + 7];
		nyE[numEny + 8].oppo = &nyE[numEny + 4];
		nyE[numEny + 8].face = &nyF[numFny + 2];

		nyF[numFny + 2].from = &nyE[numEny + 6];
		

	
		if (printar) cout << "HALLÅÅ 2: nyE[" << numEny + 5 << "] = " << nyE[numEny + 5].to - nyV << endl;
		iterE->to = nyV_[0];
		iterE->next->to = nyV_[1];
		iterE->prev->to = nyV_[2];
		if (printar) cout << "HALLÅÅ 3: nyE[" << numEny + 5 << "] = " << nyE[numEny + 5].to - nyV << endl;


		iterE->prev->face = &nyF[numFny + 1];	
		iterE->prev->prev = &nyE[numEny + 5];
		iterE->prev->next = &nyE[numEny + 4];

		iterE->next->face = &nyF[numFny + 0];
		iterE->next->prev = &nyE[numEny + 3];
		iterE->next->next = &nyE[numEny + 2];

		iterE->prev = &nyE[numEny + 1];
		iterE->next = &nyE[numEny + 0];
		

		numEny += 9;
		numFny += 3;


			////////////////////////////////////////////////
		if (printar) {
			cout << "antalV: " << numV << "\t";
			cout << "antalE: " << numE << "\t";
			cout << "antalF: " << numF << endl;
		}

		if (printar) {
			cout << "\tV[n] = {Edge *from, \tVec3 X, \tVec3 Norm}" << endl;
			for (int v_=0; v_<numVny; v_++)
				cout << "nyV[" << v_ << "] = {" << (nyV[v_].from - nyE) << ",\t" << nyV[v_].X << ",\t" << nyV[v_].Norm << endl;

			cout << endl << "\tE[n] {Vertex *fr,\tVertex *to,\tEdge *next,\tEdge *prev,\tEdge *oppo,\tFace *face}" << endl;
			for (int e_=0; e_<numEny; e_++) {
				cout << "nyE[" << e_ << "] {" << (nyE[e_].fr - nyV) << ",\t" << (nyE[e_].to - nyV) << ",\t" << (nyE[e_].next - nyE) << ",\t";
				cout << (nyE[e_].prev - nyE) << ",\t" << (nyE[e_].oppo?(nyE[e_].oppo - nyE) :-1) << ",\t" << (nyE[e_].face - nyF) << "}" << endl;
			}

			cout << endl << "\tF[n] {Edge *from,\tVec Norm}" << endl;
			for (int f_=0; f_<numFny; f_++)
				cout << "nyF[" << f_ << "] {" << nyF[f_].from - nyE << ", " << nyF[f_].Norm << "}" << endl;
		}
			/////////////////////////////////////////////////

	}
	/////////////////////////////////////////////////////////
	if (printar) cout << "Blev färdig, nu bara fixa det sista" << endl;
	delete[] E;
	delete[] F;
	delete[] V;

	E = nyE;
	F = nyF;
	V = nyV;

	numV = numVny;
	numE = numEny;
	numF = numFny;

	for (int e=0; e<numE; e++)
		E[e].fr->from = &E[e];

	
	normalizeRadius();

	cout << "numV: " << numV << endl;
	cout << "numE: " << numE << endl;
	cout << "numF: " << numF << endl;


	return true;
}


bool ObjectFN::subdivide2(int N)	// divides every edge n times. subdivide2() = subdivide(2)
{
	if (!consistsOfOnlyTriangles) {
		cout << "kan inte subdivida om inte enbart trianglar i polyeder" << endl;
		return false;
	}

	if (N < 3)
	{
		cout << "fail om N < 3. N = " << N << endl;
		return false;
	}

	bool printar = false;

	Vec _Center = getCenter();
	TYP _rad = sqrt((V[0].X - _Center) * (V[0].X - _Center));
	cout << "Center: " << _Center << endl;
	cout << "radius: " << _rad << endl;

	if (N < 2) {
		cout << "subdivide2-error, 2 <= n" << endl;
	} else if (N == 2) {
		cout << "úse funktionen subdivide(void) om vid simpel subdivision" << endl;
		return subdivide2();
	}

	int n2 = N*N;

	int numVny = numV;
	int numEny = numE;
	int numFny = numF;


	Vertex *nyV = new Vertex[2 + numF*n2/2];
	Edge *nyE = new Edge[numE * n2];
	Face *nyF = new Face[numF * n2];


	CopyVEF(nyV, nyE, nyF);

	for (int e=0; e<numE; e++)
	{
			// om det redan är fixat med denna edgen: gå vidare
		if (nyE[e].oppo - nyE >= numE)
			continue;

		Edge *_denna = &nyE[numEny];
		Edge *_motsatta = &nyE[numEny + N - 1];

		for (int n=0; n<N-2; n++)
		{
			_denna->next = _denna + 1;
			_denna->oppo = _motsatta + N - 3 - (2*n);
			_denna->prev = _denna - 1;
			_denna->face = nyE[e].face;

			_denna->fr = &nyV[numVny + n];
			_denna->fr->from = _denna;
			_denna->to = &nyV[numVny + n + 1];

			_motsatta->fr = &nyV[numVny + N - n - 2];
			_motsatta->to = &nyV[numVny + N - n - 3];
			_motsatta->fr->from = _motsatta;

			TYP frac_ = (n+1.0) / N;


			nyV[numVny + n].X = nyE[e].to->X * frac_ + nyE[e].fr->X * (1-frac_);
			nyV[numVny + n].from = _denna;


			_motsatta->next = _motsatta + 1;
			_motsatta->oppo = _denna + N - 3 - (2*n);
			_motsatta->prev = _motsatta - 1;
			_motsatta->face = nyE[e].oppo->face;

				//fixa vertices också
			_denna = _denna+1;
			_motsatta = _motsatta+1;
		}
		nyV[numVny + N-2].X = nyE[e].to->X * (1.-1./N) + nyE[e].fr->X * (1. / N);

		_motsatta->oppo 	= &nyE[e];
		_motsatta->next 	= nyE[e].oppo->next;
		_motsatta->prev 	= _motsatta - 1;
		_motsatta->fr 		= _motsatta->prev->to;
		_motsatta->fr->from = _motsatta;
		_motsatta->to 		= _motsatta->oppo->fr;
		_motsatta->face 	= nyE[e].oppo->face;
		_motsatta->next->prev = _motsatta;

		_denna->oppo 		= nyE[e].oppo;
		_denna->next 		= nyE[e].next;
		_denna->prev 		= _denna-1;
		_denna->fr 			= _denna->prev->to;
		_denna->fr->from 	= _denna;
		_denna->to 			= _denna->oppo->fr;
		_denna->face 		= nyE[e].face;
		_denna->fr->from	= _denna;
		_denna->next->prev 	= _denna;

		_motsatta->prev->oppo->prev	= &nyE[e];
		_denna->prev->oppo->prev = nyE[e].oppo;

		nyE[e].oppo->oppo 	= _denna;
		nyE[e].oppo->next 	= _denna->prev->oppo;
		nyE[e].oppo->to 	= _denna->fr;

		nyE[e].oppo 		= _motsatta;
		nyE[e].next 		= _motsatta->prev->oppo;
		nyE[e].to 			= _motsatta->fr;


		numVny += (N-1);
		numEny += 2*(N-1);
	}


	int k_;
	for (int f=0; f<numF; f++)
	{
		Vec dAB = nyF[f].from->to->X - nyF[f].from->fr->X;
		Vec dBC = nyF[f].from->prev->fr->X - nyF[f].from->to->X;

		for (int r=1; r<N; r++)
		{
			k_ = numEny + 3*r*r - 2*r;
			for (int n=0; n<2*r+1; n++)
			{
					// fixa triangel r^2 + n

				if (n&1) {
						// fixa vertex coords
					if (r != N-1)
						truncatedVertexNum(nyV + numVny, N, r, n, 1)->X = nyF[f].from->fr->X + dAB*(r+1.) + dBC*(n+1.)/2;

					nyE[k_].fr 		= truncatedVertexNum(nyV + numVny, N, r, n, 0);
					if (nyE[k_].fr)
						nyE[k_].fr->from = &nyE[k_];
					nyE[k_].to 		= truncatedVertexNum(nyV + numVny, N, r, n, 1);
					nyE[k_].next 	= truncatedEdgeNum(nyE + numEny, N, r, n, 1);
					nyE[k_].prev 	= truncatedEdgeNum(nyE + numEny, N, r, n, 2);
					nyE[k_].oppo 	= truncatedEdgeNum(nyE + numEny, N, r, n-1, 2);
					nyE[k_].face 	= &nyF[numFny + r*r - 1 + n];
					nyE[k_].face->from = &nyE[k_];
					k_++;

					nyE[k_].fr 		= truncatedVertexNum(nyV + numVny, N, r, n, 1);
					if (nyE[k_].fr)
						nyE[k_].fr->from = &nyE[k_];
					nyE[k_].to 		= truncatedVertexNum(nyV + numVny, N, r, n, 2);
					nyE[k_].next 	= truncatedEdgeNum(nyE + numEny, N, r, n, 2);
					nyE[k_].prev 	= truncatedEdgeNum(nyE + numEny, N, r, n, 0);
					nyE[k_].oppo 	= truncatedEdgeNum(nyE + numEny, N, r, n+1, 0);
					nyE[k_].face 	= &nyF[numFny + r*r - 1 + n];
					nyE[k_].face->from = &nyE[k_];
					k_++;

					nyE[k_].fr 		= truncatedVertexNum(nyV + numVny, N, r, n, 2);
					if (nyE[k_].fr)
						nyE[k_].fr->from = &nyE[k_];
					nyE[k_].to 		= truncatedVertexNum(nyV + numVny, N, r, n, 0);
					nyE[k_].next 	= truncatedEdgeNum(nyE + numEny, N, r, n, 0);
					nyE[k_].prev 	= truncatedEdgeNum(nyE + numEny, N, r, n, 1);
					nyE[k_].oppo 	= truncatedEdgeNum(nyE + numEny, N, r-1, n-1, 1);
					nyE[k_].face 	= &nyF[numFny + r*r - 1 + n];
					nyE[k_].face->from = &nyE[k_];
					k_++;


				} else {
					if (n > 0) {
						nyE[k_].fr 		= truncatedVertexNum(nyV + numVny, N, r, n, 2);
						if (nyE[k_].fr)
							nyE[k_].fr->from = &nyE[k_];
						nyE[k_].to 		= truncatedVertexNum(nyV + numVny, N, r, n, 0);
						nyE[k_].next 	= truncatedEdgeNum(nyE + numEny, N, r, n, 1);
						nyE[k_].prev 	= truncatedEdgeNum(nyE + numEny, N, r, n, 2);
						nyE[k_].oppo 	= truncatedEdgeNum(nyE + numEny, N, r, n-1, 1);
						nyE[k_].face 	= &nyF[numFny + r*r - 1 + n];
						nyE[k_].face->from = &nyE[k_];
						k_++;
					}

					if (r != N-1) {
						nyE[k_].fr 		= truncatedVertexNum(nyV + numVny, N, r, n, 0);
						if (nyE[k_].fr)
							nyE[k_].fr->from = &nyE[k_];
						nyE[k_].to 		= truncatedVertexNum(nyV + numVny, N, r, n, 1);
						nyE[k_].next 	= truncatedEdgeNum(nyE + numEny, N, r, n, 2);
						nyE[k_].prev 	= truncatedEdgeNum(nyE + numEny, N, r, n, 0);
						nyE[k_].oppo 	= truncatedEdgeNum(nyE + numEny, N, r+1, n+1, 2);
						nyE[k_].face 	= &nyF[numFny + r*r - 1 + n];
						nyE[k_].face->from = &nyE[k_];
						k_++;
					}

					if (n < 2*r)
					{
						nyE[k_].fr 		= truncatedVertexNum(nyV + numVny, N, r, n, 1);
						if (nyE[k_].fr)
							nyE[k_].fr->from = &nyE[k_];
						nyE[k_].to 		= truncatedVertexNum(nyV + numVny, N, r, n, 2);
						nyE[k_].next 	= truncatedEdgeNum(nyE + numEny, N, r, n, 0);
						nyE[k_].prev 	= truncatedEdgeNum(nyE + numEny, N, r, n, 1);
						nyE[k_].oppo 	= truncatedEdgeNum(nyE + numEny, N, r, n+1, 0);
						nyE[k_].face 	= &nyF[numFny + r*r - 1 + n];
						nyE[k_].face->from = &nyE[k_];
						k_++;
					}
				}
			}
		}

		Edge *_denna = nyF[f].from;
		Edge *_dennaNext = _denna->next;
		Edge *_motsatta = _denna->prev;
		Edge *_motsattaPrev = _motsatta->prev;

			// fixa första triangeln r=0, n=0:
		_motsatta->prev = _denna->next = truncatedEdgeNum(nyE + numEny, N, 0, 0, 1);
		_motsatta->face = _denna->face;

		_denna->next->fr = _denna->to;
		_denna->next->to = _motsatta->fr;
		_denna->next->next = _motsatta;
		_denna->next->prev = _denna;
		_denna->next->oppo = truncatedEdgeNum(nyE + numEny, N, 1, 1, 2);
		_denna->next->face = _denna->face;

			// uppdatera denna och motsatta till rad 1
		_denna = _dennaNext;
		_dennaNext = _denna->next;
		_motsatta = _motsattaPrev;
		_motsattaPrev = _motsatta->prev;

		for (int r=1; r<N-1; r++)
		{

			if (printar)	cout << "r = " << r << endl;

				// fixa vänstersidan på rad r
			_denna->next = truncatedEdgeNum(nyE + numEny, N, r, 0, 1);
			_denna->next->prev = _denna;
			_denna->next->fr = _denna->to;

			_denna->prev = truncatedEdgeNum(nyE + numEny, N, r, 0, 2);
			_denna->prev->next = _denna;
			_denna->prev->to = _denna->fr;

			_denna->prev->oppo->fr = _denna->fr;
			_denna->prev->oppo->prev->to = _denna->fr;
			_denna->face = _denna->next->face;
			if (printar)	cout << "steg 1 "<< endl;

				// fixa högersidan på rad r
			_motsatta->prev = truncatedEdgeNum(nyE + numEny, N, r, 2*r, 1);
			_motsatta->prev->next = _motsatta;
			_motsatta->prev->to = _motsatta->fr;

			_motsatta->next = truncatedEdgeNum(nyE + numEny, N, r, 2*r, 0);
			_motsatta->next->prev = _motsatta;
			_motsatta->next->fr = _motsatta->to;
			if (printar)	cout << "steg 2 "<< endl;
			
			_motsatta->next->oppo->to = _motsatta->to;
			_motsatta->next->oppo->next->fr = _motsatta->to;
			_motsatta->face = _motsatta->next->face;
			if (printar)	cout << "steg 3 "<< endl;

				// stega vidare
			_denna = _dennaNext;
			_dennaNext = _denna->next;
			_motsatta = _motsattaPrev;
			_motsattaPrev = _motsatta->prev;
		}	// end for r

			// sista raden av korrigering r= N-1,  n=0
		_denna = _dennaNext;
		_dennaNext = _denna->next;

		_denna->face = &nyF[numFny + (N-1)*(N-1) - 1];
		_denna->prev->face = _denna->face;

		_motsatta->face = _denna->face + 2*(N-1);
		_motsatta->prev->face = _motsatta->face;
		_motsatta->next = truncatedEdgeNum(nyE + numEny, N, N-1, 2*(N-1), 0);

		_denna->next = truncatedEdgeNum(nyE + numEny, N, N-1, 0, 2);
		_denna->next->prev = _denna;
		_denna->next->next = _denna->prev;
		_denna->prev->prev = _denna->next;
		_denna->next->to = _denna->next->oppo->fr = _denna->next->oppo->prev->to = _denna->prev->fr;
		_denna->next->fr = _denna->next->oppo->to = _denna->next->oppo->next->fr = _denna->to;

		_denna = _dennaNext;
		_dennaNext = _denna->next;

		for (int n=2; n<2*(N-1); n += 2)
		{
			_denna->next = truncatedEdgeNum(nyE + numEny, N, N-1, n, 2);
			_denna->prev = truncatedEdgeNum(nyE + numEny, N, N-1, n, 0);
			_denna->next->prev = _denna->prev->next = _denna;
			_denna->face = &nyF[numFny + (N-1)*(N-1) + n - 1];

			_denna->prev->to = _denna->fr;
			_denna->next->fr = _denna->next->oppo->to = _denna->next->oppo->next->fr = _denna->to;

			_denna = _dennaNext;
			_dennaNext = _denna->next;			
		}

		_denna->prev = truncatedEdgeNum(nyE + numEny, N, N-1, 2*N-2, 0);
		_denna->prev->next = _denna;
		_denna->prev->prev = _denna->next;
		_denna->prev->fr = _denna->prev->oppo->to = _denna->prev->oppo->next->fr = _denna->next->to;
		_denna->prev->to = _denna->fr;

		numVny += (N-1)*(N-2)/2;
		numEny += 3*N*N - 3*N;
		numFny += (N*N - 1);
	}	// end for f


	for (int v=0; v<numVny; v++)
	{
		nyV[v].X -= _Center;
		nyV[v].X *= _rad / sqrt(nyV[v].X*nyV[v].X);
		nyV[v].X += _Center;
	}


	for (int f=0; f<numFny; f++)
	{
		nyF[f].update();
	}

	delete[] V;
	delete[] E;
	delete[] F;

	numV = numVny;
	numE = numEny;
	numF = numFny;

	V = nyV;
	E = nyE;
	F = nyF;

	return false;
}


bool ObjectFN::makeDual()
{
	bool printar = false;

		// skapa nya vertex, edges och faces som ska användas.
	int numVny = numF;
	int numEny = numE;
	int numFny = numV;

	Vertex *nyV = new Vertex[numVny];
	Edge *nyE = new Edge[numEny];
	Face *nyF = new Face[numFny];

	for (int f=0; f<numF; f++)
	{
		if (printar) cout << "v[" <<f << "]: ";
		int i=0;
		Edge *e = F[f].from;
		Edge *endingE = e->prev;
		Vec _Pos = endingE->fr->X;
		while (e != endingE)
		{
			_Pos += e->fr->X;
			e = e->next;
			i++;
			if (i>20)
			{
				if (printar)	cout << "gick dåligt" << endl;
				return false;
			}
		}

		_Pos /= i;
		if (printar)	cout << "i: " << i << "\tP: " << nyV[f].X << endl;
		nyV[f].X = _Pos;
	}

	if (printar) 	cout << "Vert is done. " << endl;

	for (int e=0; e<numE; e++)
	{
		if (printar) 	cout << "e: " << e << endl;
		
		nyE[e].fr 	= nyV + (E[e].face - F);
		nyE[e].fr->from = &nyE[e];
		nyE[e].to 	= nyV + (E[e].oppo->face - F);
		nyE[e].next = nyE + (E[e].oppo->prev - E);
		nyE[e].prev = nyE + (E[e].next->oppo - E);
		nyE[e].oppo = nyE + (E[e].oppo - E);
		nyE[e].face = nyF + (E[e].to - V);
		nyE[e].face->from = &nyE[e];
	}

	for (int f=0; f<numFny; f++)
	{
		if (printar) 	cout << "f: " << f << endl;
		nyF[f].update();
	}


	delete[] V;
	delete[] E;
	delete[] F;
	numV = numVny;
	numE = numEny;
	numF = numFny;

	V = nyV;
	E = nyE;
	F = nyF;

	if (printar)
		cout << (updateConsistsOfOnlyTriangles()? "bara trianglar": "inte bara trianglar") << endl;
	else
		updateConsistsOfOnlyTriangles();

	return true;
}



bool ObjectFN::truncate(TYP val) 	// truncated = 0.5, rectified = 1.0;
{
	bool printar = false;

		// skapa nya vertex, edges och faces som ska användas.
	int numVny = numV;
	int numEny = numE;
	int numFny = numF;

	Vertex *nyV = new Vertex[numE];
	Edge *nyE = new Edge[2*(numV + numE + numF) - 4];
	Face *nyF = new Face[numF + numV];


	if (printar)	cout << "nyV: " << nyV << endl;
	if (printar)	cout << "nyE: " << nyE << endl;
	if (printar)	cout << "nyF: " << nyF << endl;

	if (printar)	cout << "numVny: " << numE << endl;
	if (printar)	cout << "numEny: " << 2*(numV + numE + numF) - 4 << endl;
	if (printar)	cout << "numFny: " << numF + numV << endl;

	CopyVEF(nyV, nyE, nyF);

	for (int v=0; v<numV; v++)
	{
		int _varning = 0;
		Edge *itE = nyV[v].from;

		nyE[numEny].fr = &nyV[numVny];
		nyE[numEny].fr->from = &nyE[numEny];
		nyE[numEny].to = &nyV[v];
		nyE[numEny].next = itE;
		nyE[numEny].prev = itE->prev;
		nyE[numEny].oppo = &nyE[numEny+1];
		nyE[numEny].face = itE->face;
		nyE[numEny].face->from = &nyE[numEny];
		

		nyE[numEny+1].fr = &nyV[v];
		nyE[numEny+1].to = &nyV[numVny];
		nyE[numEny+1].next = &nyE[numEny+3];
		nyE[numEny+1].oppo = &nyE[numEny];
		nyE[numEny+1].face = &nyF[numFny];
		nyE[numEny+1].face->from = &nyE[numEny+1];


		
		Edge *oldE = itE;
		
		itE = itE->prev;
		itE->next->prev = &nyE[numEny];
		itE->next = &nyE[numEny];
		itE = itE->oppo;


		do {
			numEny += 2;
			itE->fr = &nyV[numVny];
			itE->fr->X = V[v].X;

				// gammal surface
			nyE[numEny].fr = &nyV[numVny+1];
			if (nyE[numEny].fr < nyV + numE)
				nyE[numEny].fr->from = &nyE[numEny];
			//nyE[numEny].fr->from = &nyE[numEny];
			nyE[numEny].to = &nyV[numVny];
			nyE[numEny].next = itE;
			nyE[numEny].prev = itE->prev;
			nyE[numEny].oppo = &nyE[numEny+1];
			nyE[numEny].face = itE->face;
			nyE[numEny].face->from = &nyE[numEny];

				// ny surface
			nyE[numEny+1].fr = &nyV[numVny];
			nyE[numEny+1].fr->from = &nyE[numEny+1];
			nyE[numEny+1].to = &nyV[numVny+1];
			nyE[numEny+1].next = &nyE[numEny+3];
			nyE[numEny+1].prev = &nyE[numEny-1]; 
			nyE[numEny+1].prev->next = &nyE[numEny+1]; 
			nyE[numEny+1].oppo = &nyE[numEny];
			nyE[numEny+1].face = &nyF[numFny];
			nyE[numEny+1].face->from = &nyE[numEny+1];

			Edge *oldE = itE;
			itE->oppo->to = itE->fr;

			itE = itE->prev;
			itE->next->prev = &nyE[numEny];
			itE->next = &nyE[numEny];
			itE = itE->oppo;
			numVny += 1;

			if (_varning++ > 10){
				cout << "Fastnade i eternal loopness of fire death mist steeeel, yeeeeääähh " << endl;
				return false;
			}
		} while(itE != nyV[v].from);

		itE->prev->oppo->prev = itE->oppo->next->oppo;
		itE->oppo->next->oppo->next = itE->prev->oppo;
		itE->oppo->next->fr = &nyV[v];
		itE->oppo->next->oppo->to = &nyV[v];

		numEny += 2;
		numFny++;
	}

	for (int e=0; e<numE; e++)
	{
		if (&nyE[e] < nyE[e].oppo){
			Vec dX = nyE[e].to->X - nyE[e].fr->X;
			nyE[e].fr->X += dX * (val*.5);
			nyE[e].to->X -= dX * (val*.5);
		}
	}

	if (printar)	cout << "numVny: " << numVny << endl;
	if (printar)	cout << "numEny: " << numEny << endl;
	if (printar)	cout << "numFny: " << numFny << endl;



	if (printar)	cout << "old V: " << V << endl;
	if (printar)	cout << "old E: " << E << endl;
	if (printar)	cout << "old F: " << F << endl;

	delete[] V;
	delete[] E;
	delete[] F;

	V = nyV;
	E = nyE;
	F = nyF;

	for (int f=0; f<numFny; f++)
		nyF[f].update();

	numV = numVny;
	numE = numEny;
	numF = numFny;
	
	return true;
}


bool ObjectFN::rectify()
{
	bool printar = false;

	int numVny = numE/2;
	int numEny = numE;
	int numFny = numF;

	Vertex *nyV = new Vertex[numE/2];
	Edge *nyE = new Edge[numE + 2*numV + 2*numF - 4];
	Face *nyF = new Face[numV + numF];

	EmptyVEF(nyV, nyE, nyF);

	Vertex **VfrE = new Vertex*[numE];


		///////////////////////////
	for (int e=0; e<numE; e++)
		VfrE[e] = 0;

	int v = 0;
	for (int e=0; e<numE; e++)
	{
		if (!VfrE[e]) {
			VfrE[e] = &nyV[v];
			VfrE[E[e].oppo - E] = &nyV[v++];
			VfrE[e]->X = (E[e].fr->X + E[e].to->X)*.5;
		}
	}


	for (int v=0; v<numV; v++)
	{
		Edge *iterE = V[v].from;
		if (printar)	cout << "iterE = E[" << V[v].from - E << "]" << endl;

		int numEround = 0;		// Antalet vertieces runt hörnet
		Edge *nyYttreE = &nyE[iterE->oppo - E];
		Edge *nyInreE = &nyE[numEny];

		do {
			nyYttreE->fr = VfrE[iterE - E];
			nyYttreE->to = VfrE[iterE->oppo->next - E];
			nyYttreE->next = &nyE[iterE->oppo->next - E];
			nyYttreE->prev = &nyE[iterE->oppo->prev - E];
			nyYttreE->oppo = nyInreE;
			nyYttreE->face = &nyF[iterE->oppo->face - F];
			nyYttreE->face->from = nyYttreE;


			nyInreE->fr = VfrE[iterE->oppo->next - E];
			nyInreE->fr->from = nyInreE;
			nyInreE->to = VfrE[iterE - E];
			nyInreE->prev = (iterE->oppo->next == V[v].from)? (&nyE[numEny + 0]): (nyInreE + 1);nyInreE->prev->next = nyInreE;
			nyInreE->oppo = nyYttreE;
			nyInreE->face = &nyF[numF + v];
			nyInreE->face->from = nyInreE;

			iterE = iterE->oppo->next;

			nyInreE++;
			nyYttreE = &nyE[iterE->oppo - E];

			numEround++;
		} while(iterE != V[v].from);

		numEny += numEround;
		numFny += 1;
	}

	for (int f=0; f<numFny; f++)
		nyF[f].update();

	delete[] V;
	delete[] E;
	delete[] F;

	numV = numVny;
	numE = numEny;
	numF = numFny;

	V = nyV;
	E = nyE;
	F = nyF;

	delete[] VfrE;
	return true;
}


bool ObjectFN::expand(TYP val) {		// val är en radiella förändringsfaktorn, val > 1.


	bool printar = false;

	int numVny = numV;
	int numEny = numE;
	int numFny = numF;

	Vertex *nyV = new Vertex[2 + 3*numE/2 - numF - numV];
	Edge *nyE = new Edge[4*numE];
	Face *nyF = new Face[numV + numE/2 + numF];

	if (printar)		cout << "numVny: " << 2 + 3*numE/2 - numF - numV << endl;
	if (printar)		cout << "numEny: " << 4*numE << endl;
	if (printar)		cout << "numFny: " << numV + numE/2 + numF << endl;


	frV__ = nyV;
	toV__ = nyV + 2 + 3*numE/2 - numF - numV;
	frE__ = nyE;
	toE__ = nyE + 4*numE;
	frF__ = nyF;
	toF__ = nyF + numV + numE/2 + numF;




	CopyVEF(nyV, nyE, nyF);


	for (int v=0; v<numV; v++)
	{
		if (printar)		cout << endl << "\tv: " << v << endl;
		Edge *iterE = nyV[v].from;
		if (printar)		cout << "iterE =nyE[" << iterE - nyE << "]" << endl;
		check(iterE, 0);
		nyV[v].X += iterE->face->Norm*val;
		if (printar)		cout << "nyV[" << v << "].X = " << nyV[v].X << endl;
		iterE = iterE->oppo->next;
		if (printar)		cout << "iterE = " << iterE - nyE << endl;
		check(iterE, 0);
		int k= 0;
		Edge *firstEdge = &nyE[numEny];
		do {

				// befintliga sidor och kanter uppdateras
			iterE->fr = iterE->prev->to = &nyV[numVny];
			if (printar)		cout << "iterE->fr = nyE[" << iterE - nyE << "].fr = iterE->prev->to = nyE[" << iterE->prev - nyE << "].to = nyV[" << numVny << "]" << endl;
			check(iterE->fr, 1);
			
				// ny vertex skapas
			//nyV[numVny].X = iterE->fr->X;
			nyV[numVny].X = V[v].X;
			nyV[numVny].X += iterE->face->Norm*val;
			if (printar)		cout << "nyV[" << numVny << "].X = " << V[v].X << " + " << iterE->face->Norm*val << " = " << nyV[numVny].X << endl;

			nyV[numVny].from = &nyE[numEny];
			if (printar)		cout << "nyV[" << numVny << "].from = nyE[" << numEny << "]" << endl;
			check(nyV[numVny].from, 2);

			//iterE->prev->to = iterE->fr = &nyV[numVny];
			//cout << "nyE[" << iterE->prev - nyE << "].to = nyE[" << iterE - nyE << "].fr = nyV[" << numVny << "]" << endl;
			//check(iterE->fr, 3);



				// sida skapad vid kanten uppdateras
			nyE[numEny].fr = iterE->fr;
			if (printar)		cout << "nyE[" << numEny << "].fr = nyV[" << iterE->fr - nyV << "]" << endl;
			check(nyE[numEny].fr, 4);
			nyE[numEny].to = iterE->prev->oppo->fr;
			if (printar)		cout << "nyE[" << numEny << "].to = nyV[" << iterE->prev->oppo->fr  - nyV<< "]" << endl;
			check(nyE[numEny].to, 5);
			nyE[numEny].face = &nyF[numFny];
			if (printar)		cout << "nyE[" << numEny << "].face = nyF[" << numFny << "]" << endl;
			check(nyE[numEny].face, 6);

			nyE[numEny].prev = &nyE[numEny+1];
			if (printar)		cout << "nyE[" << numEny << "].prev = nyE[" << numEny+1 << "]" << endl;
			nyE[numEny].prev->next = &nyE[numEny];
			if (printar)		cout << "nyE[" << nyE[numEny].prev - nyE << "].next = nyE[" << numEny << "]" << endl;



			numVny++;
			numEny++;
			if (printar)		cout << "numXny = " << numVny << ", " << numEny << ", " << numFny << endl;
			iterE = iterE->oppo->next;
			if (printar)		cout << "iterE = nyE[" << iterE - nyE << "]" << endl;
			check(iterE, 7);
			if (++k>5){
				cout << "k blev too big" << endl;
				return false;
			}

		} while(iterE != nyV[v].from);

		nyE[numEny].fr = &nyV[v];
		if (printar)		cout << "nyE[" << numEny << "].fr = nyV[" << v << "]" << endl;
		nyE[numEny].to = iterE->prev->oppo->fr;
		if (printar)		cout << "nyE[" << numEny << "].to = nyV[" << iterE->prev->oppo->fr - nyV << "]" << endl;
		nyE[numEny].face = &nyF[numFny];
		if (printar)		cout << "nyE[" << numEny << "].face = nyF[" << &nyF[numF] - nyF << "]" << endl;

		
		
		nyE[numEny].prev = firstEdge;
		if (printar)		cout << "nyE[" << numEny << "].prev = nyE[" << firstEdge - nyE << "]" << endl;
		firstEdge->next = &nyE[numEny];
		if (printar)		cout << "nyE[" << firstEdge - nyE << "].next = nyE[" << numEny << "]" << endl;

		nyV[v].from = &nyE[numEny];
		if (printar)		cout << "nyV[" << v << "].from = nyE[" << numEny - 1 << "]" << endl;
		nyF[numFny].from = &nyE[numEny];
		if (printar)		cout << "nyF[" << numFny << "].from = nyE[" << &nyE[numEny] - nyE << "]" << endl;
		numFny++;
		numEny++;
		if (printar)		cout << "numXny = " << numVny << ", " << numEny << ", " << numFny << endl;
	}


	if (printar)		cout << "printar allt:" << endl;
	if (printar)		print(nyV, nyE, nyF, numVny, numEny, numFny);
	if (printar)		cout << "Ska uppdatera fäjsen" << endl;

	if (printar)		cout << endl << endl << "Fixa edgesarnanarnarna" << endl;
	for (int e=0; e<numE; e++)
	{
		if (&nyE[e] < nyE[e].oppo)
			continue;

		if (printar)		cout << "Skapar face mellan nyE[" << e << "] och nyE[" << nyE[e].oppo - nyE << "]" << endl;

		nyE[numEny+0].fr = nyE[e].to;
		nyE[numEny+0].to = nyE[e].fr;
		nyE[numEny+0].next = &nyE[numEny+1];
		nyE[numEny+0].prev = &nyE[numEny+3];
		nyE[numEny+0].oppo = &nyE[e];
		nyE[numEny+0].face = &nyF[numFny];

		nyE[numEny+1].fr = nyE[e].fr;
		nyE[numEny+1].to = nyE[e].oppo->to;
		nyE[numEny+1].next = &nyE[numEny+2];
		nyE[numEny+1].prev = &nyE[numEny+0];
		nyE[numEny+1].oppo = nyE[e].fr->from->prev;
		if (printar)		cout << "nyE[" << numEny + 1 << "].oppo = nyE[" << nyE[e].fr->from->prev -nyE << "]" << endl;
		nyE[numEny+1].face = &nyF[numFny];

		nyE[numEny+2].fr = nyE[e].oppo->to;
		nyE[numEny+2].to = nyE[e].oppo->fr;
		nyE[numEny+2].next = &nyE[numEny+3];
		nyE[numEny+2].prev = &nyE[numEny+1];
		nyE[numEny+2].oppo = nyE[e].oppo;
		nyE[numEny+2].face = &nyF[numFny];

		nyE[numEny+3].fr = nyE[e].oppo->fr;
		nyE[numEny+3].to = nyE[e].to;
		nyE[numEny+3].next = &nyE[numEny+0];
		nyE[numEny+3].prev = &nyE[numEny+2];
		nyE[numEny+3].oppo = nyE[e].to->from;
		if (printar)		cout << "nyE[" << numEny + 3 << "].oppo = nyE[" << nyE[e].to->from -nyE << "]" << endl;
		nyE[numEny+3].face = &nyF[numFny];


		nyE[e].oppo->to->from->oppo = &nyE[numEny+1];
		if (printar)		cout << "nyE[" << nyE[e].oppo->to->from - nyE << "].oppo = nyE[" << numEny + 1 << "]" << endl;

		nyE[e].to->from->oppo = &nyE[numEny+3];
		if (printar)		cout << "nyE[" << nyE[e].to->from - nyE << "].oppo = nyE[" << numEny + 3 << "]" << endl;

		nyE[e].oppo->oppo = &nyE[numEny+2];
		if (printar)		cout << "nyE[" << nyE[e].oppo - nyE << "].oppo = nyE[" << numEny + 2 << "]" << endl;

		nyE[e].oppo = &nyE[numEny+0];
		if (printar)		cout << "nyE[" << e << "].oppo = nyE[" << numEny + 0 << "]" << endl;
		
		nyF[numFny].from = &nyE[numEny+0];

		numEny += 4;
		numFny++;
	}


	for (int f=0; f<numFny; f++)
	{
		nyF[f].update();
	}

	if (printar)		cout << "Kom hit " << endl;

	delete[] V;
	delete[] E;
	delete[] F;

	numV = numVny;
	numE = numEny;
	numF = numFny;


	cout << "numV: " << numV << endl;
	cout << "numE: " << numE << endl;
	cout << "numF: " << numF << endl;

	V = nyV;
	E = nyE;
	F = nyF;

	return true;
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

	cout << "numF = " << numF << endl;


		// kontrollera först att inga vertexes är bundna till två fäjses med N vertices
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
			iterE = iterE->next;


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
			nyFN = new ObjDodecahedronFN(Pos, Siz*sqrt(2. + 2./sqrt(5.)), Ori);
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

