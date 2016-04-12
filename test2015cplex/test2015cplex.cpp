// test2015cplex.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


int main()
{
	using namespace Eigen;
	const int Num_of_Probs = 7;
	VectorXd x0(2), x1(2), x2(2), x3(2), x4(2), x5(2), x6(2);
	x0 << 1, -0.1;
	x1 << 2, 2;
	x2 << 1, 1;
	x3 << -1, 5;
	x4 << -0.5, -0.5;
	x5 << 0.8, 0.6;
	x6 << 3, 2;
	Ref<VectorXd> X[Num_of_Probs] = {
		x0,x1,x2,x3,x4,x5,x6 };
	ProbData probArray[Num_of_Probs] = {
		{ "CB2",1.9522245 },
		{ "CB3",2 },
		{ "DEM",-3 },
		{ "QL",7.2 },
		{ "LQ",-sqrt(2) },
		{ "Mifflin1",-1 },
		{ "Wolfe",-8 } };
	void(*FunctionPointers[Num_of_Probs])(const Ref<const VectorXd>&, double &, Ref<VectorXd>) = {
		CB2,CB3,DEM,QL,LQ,Mifflin1,Wolfe };
	OutData out[Num_of_Probs];
	for (auto r=6;r<Num_of_Probs;++r)
	{
	out[r]= LPB(probArray[r],X[r],FunctionPointers[r]);
	}
	int n = 100;
	VectorXd clq(n);
	clq.setOnes();
	clq *= -0.5;
	ProbData pr = { "Chained_LQ",-(n - 1)*sqrt(2) };
	//LPB(pr, clq, Chained_LQ);


	double y0;
	VectorXd s0(2);
	Wolfe(x6, y0, s0);
	std::cout << "x6=" << x6 << " and y0=" << y0 << " and s0=" << s0 << std::endl;


    return 0;
}

