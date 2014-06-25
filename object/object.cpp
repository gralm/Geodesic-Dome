#include "object.hpp"
#include "objCube.hpp"
#include "objTetrahedron.hpp"
#include "objDodecahedron.hpp"
#include "objTruncatedIcosahedron.hpp"
		//objTruncatedIcosahedron

#define ERR_PRINT(str_) 	cout << numOfErrs << ": " << str_ << endl; if (numOfErrs++ >= maxNumOfErrs) return false

using namespace std;

std::list<ObjectFN*> World::Objs;


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

void Edge::print(Vertex *V0, Edge *E0, Face *F0)
{
	cout << "E[" << (this - E0) << "] {fr=" << (fr-V0) << ",\tto=" << (to-V0);
	cout << ",\tnext=" << (next-E0) << ",\tprev=" << (prev-E0) << ",\toppo=" << (oppo-E0);
	cout << ",\tface=" << (face-F0) << "}" << endl;
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

	V = new Vertex[_vert];
	E = new Edge[_edge];
	F = new Face[_face];
	cout << "V: " << V << endl;
	cout << "E: " << E << endl;
	cout << "F: " << F << endl;

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

void ObjectFN::CopyVEF(Vertex *nyV, Edge *nyE, Face *nyF)
{
	for (int v=0; v<numV; v++)
	{
		nyV[v].from = nyE + (V[v].from - E);
		nyV[v].X = V[v].X;
		nyV[v].Norm = V[v].Norm;
	}

	for (int e=0; e<numE; e++)
	{
		nyE[e].fr = nyV + (E[e].fr - V);
		nyE[e].to = nyV + (E[e].to - V);
		nyE[e].next = nyE + (E[e].next - E);
		nyE[e].prev = nyE + (E[e].prev - E);
		nyE[e].oppo = nyE + (E[e].oppo - E);
		nyE[e].face = nyF + (E[e].face - F);
	}

	for (int f=0; f<numF; f++)
	{
		nyF[f].from = nyE + (F[f].from - E);
		nyF[f].Norm = F[f].Norm;
	}
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

		if (edgeLen < minEdgeLenSq) {
			minEdgeLenSq = edgeLen;
			minEdgeLenId = e;
		}
		
		if (edgeLen < minEdgeLenSq) {
			minEdgeLenSq = edgeLen;
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


	cout << "nyV: " << nyV << endl;
	cout << "nyE: " << nyE << endl;
	cout << "nyF: " << nyF << endl;

	cout << "numVny: " << numE << endl;
	cout << "numEny: " << 2*(numV + numE + numF) - 4 << endl;
	cout << "numFny: " << numF + numV << endl;

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
			nyE[numEny].fr->from = &nyE[numEny];
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

	cout << "numVny: " << numVny << endl;
	cout << "numEny: " << numEny << endl;
	cout << "numFny: " << numFny << endl;



	cout << "old V: " << V << endl;
	cout << "old E: " << E << endl;
	cout << "old F: " << F << endl;

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


ObjectFN* World::addObjectFN(int objType, const Vec &Pos, const Vec &Siz, const Mat &Ori)
{

	ObjectFN *nyFN = 0;
	switch(objType)
	{
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

