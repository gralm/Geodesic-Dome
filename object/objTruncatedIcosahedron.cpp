#include "objTruncatedIcosahedron.hpp"

ObjTruncatedIcosahedronFN::ObjTruncatedIcosahedronFN(const Vec &Pos, const Vec &Siz, const Mat &Ori): ObjectFN(60, 180, 32) {
	int f, e, v;
	TYP s5 = static_cast<TYP>(sqrt(5.0));
	TYP sin1 = sqrt((5.0+s5) / 8.);
	TYP sin2 = sqrt((5.0-s5) / 8.);
	TYP cos1 = 0.25 * (s5 - 1);
	TYP cos2 = -0.25 * (s5 + 1);



				// A
	V[0].X = Vec(1.0, 0.0, 0.0);
	V[1].X = Vec(cos1, sin1, 0.0);
	V[2].X = Vec(cos2, sin2, 0.0);
	V[3].X = Vec(cos2, -sin2, 0.0);
	V[4].X = Vec(cos1, -sin1, 0.0);


	for (v=0; v<5; v++) {
				// B
		V[v+5].X = V[v].X * 2.0 + Vec(0, 0, -cos1*2);

				// C, C'
		V[v+10].X = V[v+5].X + V[_HI(v)].X + Vec(0, 0, -cos1*2);
		V[v+15].X = V[v+5].X + V[_LO(v)].X + Vec(0, 0, -cos1*2);

				// D, D'
		V[v+20].X = V[v].X*(1.5+.5*s5) + V[_HI(v)].X*(.5*s5-.5) + Vec(0, 0, -s5);
		V[v+25].X = V[v].X*(1.5+.5*s5) + V[_LO(v)].X*(.5*s5-.5) + Vec(0, 0, -s5);
	}

	for (v=0; v<30; v++)
		V[v].X += Vec(0, 0, .5 + s5);
 

	for (v=0; v<5; v++)
	{
				// E, E'
		V[v+30].X = -V[((v+3)%5) + 25].X;
		V[v+35].X = -V[((v+2)%5) + 20].X;

				// F, F'
		V[v+40].X = -V[((v+3)%5) + 15].X;
		V[v+45].X = -V[((v+2)%5) + 10].X;

				// G
		V[v+50].X = -V[((v+2)%5) + 5].X;

				// H
		V[v+55].X = -V[(v+2)%5].X;
	}

 

	E[0].set(V, E, F,	0,	1,	1,	4,	11,	0);
	E[1].set(V, E, F,	1,	2,	2,	0,	17,	0);
	E[2].set(V, E, F,	2,	3,	3,	1,	23,	0);
	E[3].set(V, E, F,	3,	4,	4,	2,	29,	0);
	E[4].set(V, E, F,	4,	0,	0,	3,	5,	0);
	F[0].from = &E[0];

	for (f=0; f<5; f++)
	{
		E[f*6 + 5].set(V, E, F,		f,			_LO(f),		6*f+6,	6*f+10,	_LO(f),			f+1);
		E[f*6 + 6].set(V, E, F,		_LO(f),		_LO(f)+5,	6*f+7,	6*f+5,	_LO(f)*6 + 10,	f+1);
		E[f*6 + 7].set(V, E, F,		_LO(f)+5,	_LO(f)+10,	6*f+8,	6*f+6,	_LO(f)*5 + 69,	f+1);
		E[f*6 + 8].set(V, E, F,		_LO(f)+10,	f+15,		6*f+9,	6*f+7,	f*6+35,			f+1);
		E[f*6 + 9].set(V, E, F,		f+15,		f+5,		6*f+10,	6*f+8,	f*5+65,			f+1);
		E[f*6 + 10].set(V, E, F,	f+5,		f,			6*f+5,	6*f+9,	_HI(f)*6 + 6,	f+1);
		F[f+1].from = &E[f*6 + 5];
	}

	for (f=0; f<5; f++)
	{
		E[f*6 + 35].set(V, E, F,	f+15,		_LO(f)+10,	6*f+36,	6*f+40,	f*6+8,			f+6);
		E[f*6 + 36].set(V, E, F,	_LO(f)+10,	_LO(f)+20,	6*f+37,	6*f+35,	_LO(f)*5 + 68,	f+6);
		E[f*6 + 37].set(V, E, F,	_LO(f)+20,	_LO(f)+30,	6*f+38,	6*f+36,	_LO(f)*6 + 92,	f+6);
		E[f*6 + 38].set(V, E, F,	_LO(f)+30,	f+35,		6*f+39,	6*f+37, _LO(f)*5 + 122,	f+6);
		E[f*6 + 39].set(V, E, F,	f+35,		f+25,		6*f+40,	6*f+38,	f*6+94,			f+6);
		E[f*6 + 40].set(V, E, F,	f+25,		f+15,		6*f+35,	6*f+39,	f*5+66,			f+6);
		F[f+6].from = &E[f*6 + 35];
	}


	for (f=0; f<5; f++)
	{
		E[f*5 + 65].set(V, E, F, 	f+5,	f+15,	f*5+66,		f*5+69, 	f*6+9, 			f+11);
		E[f*5 + 66].set(V, E, F, 	f+15,	f+25,	f*5+67,		f*5+65, 	f*6+40, 		f+11);
		E[f*5 + 67].set(V, E, F, 	f+25,	f+20,	f*5+68,		f*5+66, 	f*6+93, 		f+11);
		E[f*5 + 68].set(V, E, F, 	f+20,	f+10,	f*5+69,		f*5+67, 	_HI(f)*6 + 36, 	f+11);
		E[f*5 + 69].set(V, E, F, 	f+10,	f+5,	f*5+65,		f*5+68, 	_HI(f)*6 + 7, 	f+11);
		F[f+11].from = &E[f*6 + 65];
	}


	for (f=0; f<5; f++)
	{
		E[f*6 + 90].set(V, E, F,	f+45,	f+40,	6*f+91,		6*f+95,	f*6+148,		f+16);
		E[f*6 + 91].set(V, E, F,	f+40,	f+30,	6*f+92,		6*f+90,	f*5+123,		f+16);
		E[f*6 + 92].set(V, E, F,	f+30,	f+20,	6*f+93,		6*f+91,	_HI(f)*6+37,	f+16);
		E[f*6 + 93].set(V, E, F,	f+20,	f+25,	6*f+94,		6*f+92,	f*5+67,			f+16);
		E[f*6 + 94].set(V, E, F,	f+25,	f+35,	6*f+95,		6*f+93,	f*6+39,			f+16);
		E[f*6 + 95].set(V, E, F,	f+35,	f+45,	6*f+90,		6*f+94,	_LO(f)*5+121,	f+16);
		F[f+16].from = &E[f*6 + 90];
	}

	for (f=0; f<5; f++)
	{
		E[f*5 + 120].set(V, E, F, 	_HI(f)+50,	_HI(f)+45,	f*5 + 121,	f*5 + 124,	_HI(f)*6 + 149,	f+21);
		E[f*5 + 121].set(V, E, F, 	_HI(f)+45,	_HI(f)+35,	f*5 + 122,	f*5 + 120,	_HI(f)*6 + 95,	f+21);
		E[f*5 + 122].set(V, E, F, 	_HI(f)+35,	f+30,		f*5 + 123,	f*5 + 121,	_HI(f)*6 + 38,	f+21);
		E[f*5 + 123].set(V, E, F, 	f+30,		f+40,		f*5 + 124,	f*5 + 122,	f*6 + 91,		f+21);
		E[f*5 + 124].set(V, E, F, 	f+40,		_HI(f)+50,	f*5 + 120,	f*5 + 123,	f*6 + 147,		f+21);
		F[f+21].from = &E[f*5 + 120];
	}

	for (f=0; f<5; f++)
	{
		E[f*6 + 145].set(V, E, F,	f+55,			_HI(f) + 55,	6*f+146,	6*f+150,	179-_LO(f),			f+26);
		E[f*6 + 146].set(V, E, F,	_HI(f) + 55,	_HI(f) + 50,	6*f+147,	6*f+145,	_HI(f)*6 + 150,		f+26);
		E[f*6 + 147].set(V, E, F,	_HI(f) + 50,	f+40,			6*f+148,	6*f+146,	f*5 + 124,			f+26);
		E[f*6 + 148].set(V, E, F,	f+40,			f+45,			6*f+149,	6*f+147,	f*6+90,				f+26);
		E[f*6 + 149].set(V, E, F,	f+45,			f+50,			6*f+150,	6*f+148,	_LO(f)*5 + 120,		f+26);
		E[f*6 + 150].set(V, E, F,	f+50,			f+55,			6*f+145,	6*f+149,	_LO(f)*6 + 146,		f+26);
		F[f+26].from = &E[f*6 + 145];
	}

	E[175].set(V, E, F,	56,	55,	176,	179,	145,	31);
	E[176].set(V, E, F,	55,	59,	177,	175,	169,	31);
	E[177].set(V, E, F,	59,	58,	178,	176,	163,	31);
	E[178].set(V, E, F,	58,	57,	179,	177,	157,	31);
	E[179].set(V, E, F,	57,	56,	175,	178,	151,	31);
	F[31].from = &E[175];

	for (int v=0; v<numV; v++)
		V[v].X = Vec(V[v].X.x*Siz.x, V[v].X.y*Siz.y, V[v].X.z*Siz.z) * Ori + Pos;


	for (f=0; f<numF; f++)
		F[f].update();

	consistsOfOnlyTriangles = false;
}