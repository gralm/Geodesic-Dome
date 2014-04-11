#include <cmath>
#include <iostream>

using namespace std;

//#define V(x)			(x%5)

class Vec {

	public:
	double x, y, z;

	Vec()
	{
		x = y = z = 0;
	}

	Vec(double x_, double y_, double z_)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	double operator*(Vec X)
	{
		return x*X.x + y*X.y + z*X.z;
	}

	double abs()
	{
		return sqrt(x*x + y*y + z*z);
	}

	Vec norm()
	{
		double inv_ = abs();

		return Vec(x/inv_, y/inv_, z/inv_);
	}

	Vec print()
	{
		cout << "[" << x << ", " << y << ", " << z << "]" << endl;
	}

	Vec print(const char *str,  int i)
	{
		cout << str << "[" << i << "]:\t";
		print();
	}

	Vec operator&(Vec X)
	{
		return Vec(y*X.z - z*X.y, z*X.x - x*X.z, x*X.y - y*X.x);
	}

	Vec operator-(Vec X)
	{
		return Vec(x-X.x, y-X.y, z-X.z);
	}

	Vec operator+(Vec X)
	{
		return Vec(x+X.x, y+X.y, z+X.z);
	}

	Vec operator*(double a)
	{
		return Vec(x*a, y*a, z*a);
	}

	Vec operator/(double a)
	{
		return Vec(x/a, y/a, z/a);
	}
};


Vec nextPoint(Vec P1, Vec P2, Vec P3)
{
	Vec P4 = P3 + P1 - P2;
	double a = (P3-P2) * (P2-P1);
	a /= (P3-P2).abs() * (P2-P1).abs();
	P4 = P4 + (P3 - P2)*(2*a);
	return P4;
}

int main()
{
	double f = sqrt(5);
	Vec A[5];
	Vec B[5];
	Vec C[5], Cp[5];
	Vec D[5], Dp[5];
	Vec E[5], Ep[5];
	Vec F[5], Fp[5];
	Vec N[25];

	Vec M(0, 0, .5 + f);



	for (int i=0; i<5; i++) {
		A[i] = Vec(cos(2*M_PI*i/5), sin(2*M_PI*i/5), 0);
		A[i].print("A", i);
		cout << "A-M: " << (A[i]-M).abs() << endl;
	}

	double c = A[1].x;
	double s = A[1].y;
	double h = 2*c;
	double l = sqrt((5 - f)/2);

	cout << "c: " << c << ", s: " << s << ", h: " << h << ", l: " << l << endl;
	cout << "hej = " << (c*c - s*s + h*h) / (l*l) << endl;

	cout << endl;
	for (int i=0; i<5; i++)
	{
		B[i] = A[i]*2 + Vec(0, 0, (f-1)/2);
		B[i].print("B", i);
		(B[i]-M).print("(B-M)", i);
		cout << "B-M: " << (B[i]-M).abs() << endl;
	}

	cout << endl;
	for (int i=0; i<5; i++)
	{
		C[i] = A[i]*2 + A[(i==4? 0: (i+1))] + Vec(0, 0, f-1);
		C[i].print("C", i);
		Cp[i] = A[i]*2 + A[(i==0? 4: (i-1))] + Vec(0, 0, f-1);
		Cp[i].print("C'", i);
		(C[i]-M).print("(C-M)", i);
		cout << "C-M: " << (C[i]-M).abs() << endl;
	}

	cout << endl;
	for (int i=0; i<5; i++)
	{
		D[i] = nextPoint(Cp[i], B[i], C[i]);
		D[i].print("D", i);
		Dp[i] = nextPoint(C[i], B[i], Cp[i]);
		Dp[i].print("D'", i);
		cout << "dist: " << (D[i] - Dp[i]).abs() << endl;
		(D[i]-M).print("(D-M)", i);
		cout << "D-M: " << (D[i]-M).abs() << endl;
	}

	cout << endl;
	for (int i=0; i<5; i++)
	{
		E[i] = nextPoint(Cp[i==4? 0: i+1], C[i], D[i]);
		E[i].print("E", i);
		//Ep[i] = nextPoint(C[i], Cp[i==4? 0: i+1], Dp[i==4? 0: i+1]);
		Ep[i] = nextPoint(C[i==0? 4: i-1], Cp[i], Dp[i]);
		Ep[i].print("E'", i);
		cout << "dist: " << (E[i] - Ep[i]).abs() << "\tdist2: " << (Ep[i] - Dp[i]).abs() << endl;
		(E[i]-M).print("(E-M)", i);
		cout << "E-M: " << (E[i]-M).abs() << endl;
	}


	cout << endl;
	for (int i=0; i<5; i++)
	{
		F[i] = nextPoint(Dp[i], D[i], E[i]);
		F[i].print("F", i);
		Fp[i] = nextPoint(D[i], Dp[i], Ep[i]);
		Fp[i].print("F'", i);
		cout << "dist: " << (F[i] - Fp[i]).abs() << "\tdist2: " << (F[i] - E[i]).abs() << "\tdist3: " << (Fp[i] - Ep[i]).abs() << endl;
		cout << "F-M: " << (F[i]-M).abs() << endl;
	}




	for (int i=0; i<4; i++)
	{
	//	( A[i] * (3 + f)/2 + A[i+1] * (f-1) / 2).print("D''", i);
	}

	for (int i=0; i<5; i++) {

		//N[10+i] = (Cp[i==4?0:i+1] - C[i])&(D[i] - C[i])).norm();//.print("Normal", 13);
	}

	(((A[1] - A[0]) & (A[4] - A[0])).norm()).print();

	((B[0] - A[0]) & (A[1] - A[0])).norm().print();
	cout << endl;
	for (int i=0; i<5; i++) {
		//((B[i] - A[i]) & (A[i==4?0:(i+1)] - A[i])).norm();
		N[i] = ((B[i] - A[i]) & (A[i==4?0:(i+1)] - A[i])).norm();
		N[i].print("N", i);
	}



	cout << endl;
	for (int i=0; i<5; i++) {
		//((B[i] - A[i]) & (A[i==4?0:(i+1)] - A[i])).norm();
		N[i] = ((B[i] - A[i]) & (A[i==4?0:(i+1)] - A[i])).norm();
		N[i].print();
	}

	cout << "N0 * N1 = " << N[0]*N[1] << endl;
	cout << "theta 2 = " << acos(N[0]*N[1])*180/M_PI << endl;

	double cost = ((B[0] - A[0]) & (A[1] - A[0])).norm() * Vec(0,0,1);
	cout << "N1 * z = " << cost << endl;
	cout << "theta 1 = " << acos(cost)*180/M_PI << endl;

	cout << (1-f)/sqrt(12);


	return 0;
}
