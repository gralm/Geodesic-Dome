#include "object.hpp"
//#include <iostream>
using namespace std;



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





bool ObjectFN::kleetope()
{
	return this->subdivide1();
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



bool ObjectFN::subdivide1(int n, TYP height)	// makes pyramids of surfaces with n edges
{
	list<int> fList;

	for (int f=0; f<numF; f++)
	{
		if (F[f].countEdges() == n)
			fList.push_back(f);
	}
	int numFaces = fList.size();

	//cout << "antal killar inuti: " << numFaces << endl;

	if (numFaces <= 0)
		return false;

	Vertex *nyV = new Vertex[numV + numFaces];
	Edge *nyE = new Edge[numE + 2*n*numFaces];
	Face *nyF = new Face[numF + (n-1)*numFaces];
	//cout << "[" << numV + numFaces << ", " << numE + 2*n*numFaces << ", " << numF + (n-1)*numFaces << "]" << endl;


	CopyVEF(nyV, nyE, nyF);

	//cout << "dessa killar is in fList" << endl;
	for (list<int>::iterator itI = fList.begin(); itI != fList.end(); itI++)
	{
	//	cout << "face num: " << *itI << endl;

		Edge *nextIterE, *iterE = nyF[*itI].from;

		nyV[numV].X = nyF[*itI].getCenter() + nyF[*itI].Norm * height;
		nyV[numV].from = &nyE[numE + 1];

		for (int i=0; i<n; i++)
		{
			nextIterE = iterE->next;

			if (i > 0) {
				iterE->face = &nyF[numF + i - 1];
				nyF[numF + i - 1].from = iterE;
			}

			nyE[numE + 2*i + 0].next = &nyE[numE + 2*i + 1];
			nyE[numE + 2*i + 0].prev = iterE;
			nyE[numE + 2*i + 0].oppo = (i==n-1? &nyE[numE + 1]: &nyE[numE + 2*i + 3]);
			nyE[numE + 2*i + 0].fr = iterE->to;
			nyE[numE + 2*i + 0].to = &nyV[numV];
			nyE[numE + 2*i + 0].face = iterE->face;

			nyE[numE + 2*i + 1].next = iterE;
			nyE[numE + 2*i + 1].prev = &nyE[numE + 2*i];
			nyE[numE + 2*i + 1].oppo = (i==0? &nyE[numE + 2*n - 2]: &nyE[numE + 2*i - 2]);
			nyE[numE + 2*i + 1].fr = &nyV[numV];
			nyE[numE + 2*i + 1].to = iterE->fr;
			nyE[numE + 2*i + 1].face = iterE->face;

			iterE->next = &nyE[numE + 2*i + 0];
			iterE->prev = &nyE[numE + 2*i + 1];

			iterE->face->update();

			iterE = nextIterE;
		}

		numV++;
		numE += 2*n;
		numF += (n-1);
		//cout << "[" << numV << ", " << numE << ", " << numF << "]" << endl;
	}
	cout << endl;



	delete[] V;
	delete[] E;
	delete[] F;

	V = nyV;
	E = nyE;
	F = nyF;

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




