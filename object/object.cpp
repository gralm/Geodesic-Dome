#include "object.hpp"

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
}


void Vertex::print(const Vertex *zero)
{
	//cout << "V[" << v << "] = {" << (V[v].from - V) << ",\t" << V[v].X << ",\t" << V[v].Norm;
}

void Face::print(const Face *zero)
{

}

void Edge::print(const Edge *zero)
{

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


ObjCubeFN::ObjCubeFN(const Vec &Pos, const Vec &Siz, const Mat &Ori) :ObjectFN(8, 24, 6)
{
		cout << "Ska skapa en Kub" << endl;
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
	cout << "Skapade en kub" << endl;

	for (int i=0; i<8; i++)
		cout << V[i].X << endl;
	cout << "antalet killar \n";
	print();
	subdivide1();
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

	return true;
}



	// http://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface. 
	// Catmull Clark Subdivision av ytor
bool ObjectFN::subdivide2()
{
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

	//print();

	
	normalizeRadius();

	

	//print();

	return true;
}

/*
struct Vertex {
	Edge *from;
	Vec X;
	Vec Norm;
};

struct Face {
	Edge *from;		// Första edgen
	Vec Norm;		// Face Normal
};

struct Edge {
	Vertex *fr;
	Vertex *to;
	Edge *next;
	Edge *prev;
	Edge *oppo;
	Face *face;
*/




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

	for (int i=0; i<20; i++)
		V[i].X *= .5;


	for (int i=0; i<numV; i++)
		V[i].X = Vec(V[i].X.x*Siz.x, V[i].X.y*Siz.y, V[i].X.z*Siz.z) * Ori + Pos;

	/*static int hej = 0;
	if (hej++ == 1)
		subdivide1();*/
}





//std::vector<ObjectFN*> Objs;
//public:
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
			return 0;
	}

	if (nyFN != 0)
		Objs.push_back(nyFN);

	return nyFN;
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

