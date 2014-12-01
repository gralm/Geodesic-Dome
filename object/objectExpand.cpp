#include "object.hpp"

using namespace std;


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


	if (printar)		cout << "numV: " << numV << endl;
	if (printar)		cout << "numE: " << numE << endl;
	if (printar)		cout << "numF: " << numF << endl;

	V = nyV;
	E = nyE;
	F = nyF;

	return true;
}

bool ObjectFN::chamfer(TYP height)		// samtliga vertices måste vara connectade till 3(6) edges
{
	bool printar = false;

	for (int v=0; v<numV; v++)
	{
		if (V[v].countEdges() != 3)
		{
			cout << "Not all Vertices connected to three edges" << endl;
			return false;
		}
	}

	Vertex *nyV = new Vertex[numV*4];
	Edge *nyE = new Edge[numE*2 + 6*numV];
	Face *nyF = new Face[numF + numE/2];

	CopyVEF(nyV, nyE, nyF);
	

	for (int v=0; v<numV; v++)
	{
		Edge *iterE = nyV[v].from;


		for (int j=0; j<3; j++)	{
			if (printar)		cout << v << " : " << j << endl;
			if (printar)		cout << "E[3].to = " << nyE[3].to - nyV << "\tE[2].to = " << nyE[2].to - nyV << endl;
			if (printar)		cout << "iterE = nyE[" << iterE - nyE << "]" << endl;
			Edge *nextIterE = iterE->oppo;
			if (iterE->oppo >= nyE + numE) {
				if (printar)		cout << "\tnextIterE = " << iterE->oppo-nyE << endl;
				if (printar)		cout << "\tnextIterE = " << iterE->oppo->prev-nyE << endl;
				if (printar)		cout << "\tnextIterE = " << iterE->oppo->prev->prev-nyE << endl;
				if (printar)		cout << "\tnextIterE = " << iterE->oppo->prev->prev->prev-nyE << endl;
				if (printar)		cout << "\tnextIterE = " << iterE->oppo->prev->prev->prev->oppo-nyE << endl;
				nextIterE = iterE->oppo->prev->prev->prev->oppo;
			}
			if (printar)		cout << "nextIterE = nyE[" << nextIterE - nyE << "]" << endl;

				// första vertices 
			if (printar)		cout << "\tfirsta verticen" << endl;
			if (printar)		cout << "nyE[" << iterE - nyE << "].fr = nyV[" << numV + 3*v + j << "]" << endl;
			iterE->fr = nyV + numV + 3*v + j;
			if (printar)		cout << "nyV[" << numV + 3*v + j << "].from = nyE[" << 2*numE + 6*v + 2*j << "]" << endl;
			nyV[numV + 3*v + j].from = nyE + 2*numE + 6*v + 2*j;
			nyV[numV + 3*v + j].X = nyV[v].X;
			if (printar)		cout << "nyE[" << iterE->prev - nyE << "].to = nyV[" << iterE->fr - nyV << "]" << endl;
			iterE->prev->to = iterE->fr;
			//cout << "E[3].to = " << nyE[3].to - nyV << "\tE[2].to = " << nyE[2].to - nyV << endl;

				// första edgen;
			if (printar)		cout << "\tfirsta edgen" << endl;
			if (printar)		cout << "nyE[" << iterE - nyE << "].oppo = nyE[" << iterE + numE - nyE << "]" << endl;
			iterE->oppo = iterE + numE;
			if (printar)		cout << "nyE[" << iterE->oppo - nyE << "].fr = nyV[" << iterE->to - nyV << "]" << endl;
			iterE->oppo->fr = iterE->to;
			if (printar)		cout << "nyE[" << iterE->oppo - nyE << "].to = nyV[" << iterE->fr - nyV << "]" << endl;
			iterE->oppo->to = iterE->fr;
			if (printar)		cout << "nyE[" << iterE->oppo - nyE << "].next = nyE[" << 2*numE + 6*v + 2*j << "]" << endl;
			iterE->oppo->next = nyE + 2*numE + 6*v + 2*j;
			if (printar)		cout << "nyE[" << iterE->oppo - nyE << "].oppo = nyE[" << iterE - nyE << "]" << endl;
			iterE->oppo->oppo = iterE;
			//cout << "E[3].to = " << nyE[3].to - nyV << "\tE[2].to = " << nyE[2].to - nyV << endl;
				// andra edgen;
			if (printar)		cout << "\tandra edgen" << endl;
			if (printar)		cout << "nyE[" << iterE->oppo->next - nyE << "].fr = nyV[" << iterE->fr - nyV << "]" << endl;
			iterE->oppo->next->fr = iterE->fr;
			if (printar)		cout << "nyE[" << iterE->oppo->next - nyE << "].to = nyV[" << v << "]" << endl;
			iterE->oppo->next->to = &nyV[v];
			if (printar)		cout << "nyE[" << iterE->oppo->next - nyE << "].next = nyE[" << iterE->oppo->next + 1 - nyE << "]" << endl;
			iterE->oppo->next->next = iterE->oppo->next + 1;
			if (printar)		cout << "nyE[" << iterE->oppo->next - nyE << "].prev = nyE[" << iterE->oppo - nyE << "]" << endl;
			iterE->oppo->next->prev = iterE->oppo;
			if (printar)		cout << "nyE[" << iterE->oppo->next - nyE << "].oppo = nyE[" << iterE->oppo->next + (j==0? 5: -1) - nyE << "]" << endl;
			iterE->oppo->next->oppo = iterE->oppo->next + (j==0? 5: -1);
			//cout << "E[3].to = " << nyE[3].to - nyV << "\tE[2].to = " << nyE[2].to - nyV << endl;

				// tredje edgen
			if (printar)		cout << "\ttredje edgen" << endl;
			if (printar)		cout << "nyE[" << iterE->oppo->next->next - nyE << "].fr = nyV[" << v << "]" << endl;
			iterE->oppo->next->next->fr = &nyV[v];
			if (printar)		cout << "nyE[" << iterE->oppo->next->next - nyE << "].to = nyV[" << numV + 3*v + (j==2? 0: j+1) << "]" << endl;
			iterE->oppo->next->next->to = nyV + numV + 3*v + (j==2? 0: j+1);
			if (printar)		cout << "nyE[" << iterE->oppo->next->next - nyE << "].next = nyE[" << nextIterE + numE - nyE << "]" << endl;
			iterE->oppo->next->next->next = nextIterE + numE;

			if (printar)		cout << "*nyE[" << (nextIterE + numE) - nyE << "].prev = nyE[" << iterE->oppo->next->next - nyE << "]" << endl;
			(nextIterE + numE)->prev = iterE->oppo->next->next;

			nyV[v].from = iterE->oppo->next->next;
			iterE->oppo->next->next->next->fr = iterE->oppo->next->next->to;
			//cout << "E[3].to = " << nyE[3].to - nyV << "\tE[2].to = " << nyE[2].to - nyV << endl;			 
			

			if (printar)		cout << "nyE[" << iterE->oppo->next->next - nyE << "].prev = nyE[" << iterE->oppo->next - nyE << "]" << endl;
			iterE->oppo->next->next->prev = iterE->oppo->next;
			if (printar)		cout << "nyE[" << iterE->oppo->next->next - nyE << "].oppo = nyE[" << iterE->oppo->next->next + (j==2? -5: 1) - nyE << "]" << endl;
			iterE->oppo->next->next->oppo = iterE->oppo->next->next + (j==2? -5: 1);
			//cout << "E[3].to = " << nyE[3].to - nyV << "\tE[2].to = " << nyE[2].to - nyV << endl;


			nextIterE->oppo = nextIterE + numE;
			//nextIterE->to = nyV + numV + 3*v + 1;

			iterE = nextIterE->next;
		}

	}

		
		// fixa faces också
	for (int e=numE; e<2*numE + 6*numV; e++)
		nyE[e].face = 0;

		// expandera ursprungliga ytor
	for (int f=0; f<numF; f++)
	{
		if (printar)		cout << "fixa height on vertices kring face[" << f << "]" << endl;
		nyF[f].moveUp(height);
	}

		// fixa alla Edges->face och alla face->norm
	for (int e=0; e<numE; e++)
	{
		if (nyE[e].oppo->face != 0)
			continue;

		Edge *iterE = nyE[e].oppo;
		int numOfE = 0;
		do {
			numOfE++;
			iterE->face = nyF + numF;
			if (printar)		cout << "nyE[" << iterE - nyE << "].face = nyF[" << numF << "]" << endl;
			iterE = iterE->next;
		} while(iterE != nyE[e].oppo || numOfE > 20);

		nyF[numF].Norm = (iterE->to->X - iterE->fr->X) & (iterE->prev->prev->fr->X - iterE->fr->X);
		nyF[numF].Norm.norm();
		nyF[numF].from = iterE;
		numF++;
	}

		// fixa alla nya verts positioner:
		// Metod: http://mathworld.wolfram.com/Plane-PlaneIntersection.html
	Vec _X[3], _N[3], _P, _tf;
	for (int v=0; v<numV; v++)
	{

		_X[0] = nyV[v].from->to->X;
		_X[1] = nyV[v].from->oppo->next->to->X;
		_X[2] = nyV[v].from->prev->fr->X;
		_N[0] = nyV[v].from->oppo->face->Norm;
		_N[1] = nyV[v].from->prev->oppo->face->Norm;
		_N[2] = nyV[v].from->face->Norm;
		if (printar) {
			cout << "\t\tV[" << v << "]" << endl;
			cout << "X: ";
			cout << "\t" << _X[0] << endl;
			cout << "\t" << _X[1] << endl;
			cout << "\t" << _X[2] << endl;
			cout << "N: ";
			cout << "\t" << _N[0] << endl;
			cout << "\t" << _N[1] << endl;
			cout << "\t" << _N[2] << endl;
			cout << endl;
		}
		
		//_tf
		_P = (_N[1] & _N[2]) * (_N[0]*_X[0]);
		_P += (_N[2] & _N[0]) * (_N[1]*_X[1]);
		_P += (_N[0] & _N[1]) * (_N[2]*_X[2]);
		_P /= _N[0] * (_N[1] & _N[2]);
		if (printar) {
			cout << "P: " << _P << endl;

			cout << "0 0 0 = " << ((_P - _X[0]) * _N[0]);
			cout << ", " << ((_P - _X[1]) * _N[1]);
			cout << ", " << ((_P - _X[2]) * _N[2]) << endl;
		}
		nyV[v].X = _P;
	}

	delete[] V;
	delete[] E;
	delete[] F;

	numV = numV*4;
	numE = numE*2 + 3*numV/2;
	
	if (printar)		cout << "numV = " << numV << endl;
	if (printar)		cout << "numE = " << numE << endl;
	if (printar)		cout << "numF = " << numF << endl;

	V = nyV;
	E = nyE;
	F = nyF;
}