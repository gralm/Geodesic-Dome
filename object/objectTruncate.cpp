#include "object.hpp"
//#include <iostream>
using namespace std;



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
