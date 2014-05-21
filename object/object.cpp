#include "object.hpp"
#include "objCube.hpp"
#include "objTetrahedron.hpp"
#include "objDodecahedron.hpp"
#include "objTruncatedIcosahedron.hpp"
		//objTruncatedIcosahedron

#define ERR_PRINT(str_) 	cout << str_ << endl; if (numOfErrs++ >= maxNumOfErrs) return false

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
	_V[_fr].from = this;
}

void Face::update()
{
	Norm = (from->to->X - from->fr->X) & (from->next->to->X - from->next->fr->X);
	Norm.norm();
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




bool ObjectFN::test() const
{
	int maxNumOfErrs = 30;
	int numOfErrs = 0;
	int maxNumEdgesPerFace = 10;

	/*for (int e=0; e<numE; e++)
		cout << "E[" << e << "].length = " << sqrt((E[e].to->X - E[e].fr->X) * (E[e].to->X - E[e].fr->X)) << endl;

	for (int v=0; v<numV; v++)
		cout << "V[" << v << "].length = " << sqrt(V[v].X*V[v].X) << endl;*/


	for (int v=0; v<numV; v++)
	{
		if (V[v].from < E || V[v].from >= E + numE)
			cout << "V[" << v << "].from outside legal interval" << endl;
		else if (V[v].from->fr != &V[v])
			cout << "V[" << v << "].from->fr != &V[" << v << "]" << endl;
	}
	for (int e=0; e<numE; e++)
	{
		if (E[e].fr < V || E[e].fr >= V+numV){
			ERR_PRINT("E[" << e << "].fr outside legal interval");
		}
		
		if (E[e].to < V || E[e].to >= V+numV){
			ERR_PRINT("E[" << e << "].to outside legal interval");
			//cout <<  << endl;
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


ObjectFN* World::addObjectFN(int objType, const Vec &Pos, const Vec &Siz, const Mat &Ori)
{

	ObjectFN *nyFN = 0;
	switch(objType)
	{
		case OBJ_CUBE:{
			nyFN = new ObjCubeFN(Pos, Siz, Ori);
			break;
		}
		case OBJ_TRUNCATED_ICOSAHEDRON:{
			nyFN = new ObjTruncatedIcosahedronFN(Pos, Siz, Ori);
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

