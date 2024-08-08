#include <iostream>
#include <fstream>
#include <map>
#include <array>
#include <complex>
#include "operator.h"

using namespace std;
using namespace std::complex_literals;

int main(){
	const int D = 5;
	double sgn;
	map<pair<int, int>, top<D>> op_table;
	map<pair<int, int>, Eigen::Matrix<complex<double>, D, D>> OKQ_table;
	map<array<int, 4>, double> V;
	Eigen::Matrix<complex<double>, D, D> Ham = Eigen::Matrix<complex<double>, D, D>::Zero();

	for (int K=1; K<5; ++K)
		for (int Q=-K; Q<=K; ++Q)
			op_table[make_pair(K, Q)] = top<D>(K, Q);
	
	for (int K=1; K<5; ++K){
		OKQ_table[make_pair(K, 0)] = op_table[make_pair(K, 0)].get_matrix(); 
		for (int Q=1; Q<=K; ++Q){
			sgn = Q%2==0 ? 1.0d : -1.0d;
			OKQ_table[make_pair(K, Q)] = pow(0.5d, 0.5d)*( sgn*top<D>(K, Q).get_matrix() + top<D>(K,-Q).get_matrix() );
			OKQ_table[make_pair(K, -Q)] = 1.0i*pow(0.5d, 0.5d)*( top<D>(K, Q).get_matrix() - sgn*top<D>(K,-Q).get_matrix() );
		}
	}

	for (int K1=1; K1<5; ++K1)
		for (int K2=1; K2<5; ++K2)
			for (int Q1=-K1; Q1<=K1; ++Q1)
				for (int Q2=-K2; Q2<=K2; ++Q2)
					V[array<int, 4>{K1,K2,Q1,Q2}] = 0.0d;

	// dipole-dipole	
	V[array<int, 4>{1,1,-1,-1}] = 1.62d;
	V[array<int, 4>{1,1,0,0}]   = 4.17d;
	V[array<int, 4>{1,1,1,-1}]  = 1.26d;
	V[array<int, 4>{1,1,-1,1}]  = 1.26d;
	V[array<int, 4>{1,1,1,1}]   = 1.62d;

	// quadrupole-quadrupole
	V[array<int, 4>{2,2,-2,-2}] = -0.41d;
	V[array<int, 4>{2,2,-1,-1}] = -0.79d;
	V[array<int, 4>{2,2,0,-2}]  =  0.16d;
	V[array<int, 4>{2,2,-2,0}]  =  0.16d;
	V[array<int, 4>{2,2,0,0}]   =  1.32d;
	V[array<int, 4>{2,2,1,-1}]  = -0.23d;
	V[array<int, 4>{2,2,-1,1}]  = -0.23d;
	V[array<int, 4>{2,2,1,1}]   = -0.79d;
	V[array<int, 4>{2,2,2,2}]   = -0.58d;

	// octupole-octupole
	V[array<int, 4>{3,3,-3,-3}] =  1.16d;
	V[array<int, 4>{3,3,-2,-2}] = -1.49d;
	V[array<int, 4>{3,3,-1,-3}] = -0.14d;
	V[array<int, 4>{3,3,-3,-1}] = -0.14d;
	V[array<int, 4>{3,3,-1,-1}] =  0.80d;
	V[array<int, 4>{3,3,0,-2}]  = -0.79d;
	V[array<int, 4>{3,3,-2,0}]  = -0.79d;
	V[array<int, 4>{3,3,0,0}]   =  2.35d;
	V[array<int, 4>{3,3,1,-3}]  = -0.29d;
	V[array<int, 4>{3,3,-3,1}]  = -0.29d;
	V[array<int, 4>{3,3,1,-1}]  = -0.98d;
	V[array<int, 4>{3,3,-1,1}]  = -0.98d;
	V[array<int, 4>{3,3,1,1}]   =  0.80d;
	V[array<int, 4>{3,3,2,2}]   = -1.89d;
	V[array<int, 4>{3,3,3,-1}]  =  0.29d;
	V[array<int, 4>{3,3,-1,3}]  =  0.29d;
	V[array<int, 4>{3,3,3,1}]   =  0.14d;
	V[array<int, 4>{3,3,1,3}]   =  0.14d;
	V[array<int, 4>{3,3,3,3}]   =  1.16d;

	// dipole-octupole
	V[array<int, 4>{1,3,-1,-1}] =  1.97d;
	V[array<int, 4>{1,3,-1, 1}] = -0.20d;
	V[array<int, 4>{1,3, 1,-1}] = -0.20d;
	V[array<int, 4>{1,3,-1, 3}] = -0.89d;
	V[array<int, 4>{1,3, 0,-2}] = -0.97d;
	V[array<int, 4>{1,3, 0, 0}] =  2.38d;
	V[array<int, 4>{1,3, 1,-3}] =  0.89d;
	V[array<int, 4>{1,3, 1, 1}] =  1.97d;
	V[array<int, 4>{3,1,-1,-1}] =  1.97d;
	V[array<int, 4>{3,1,-1, 1}] = -0.20d;
	V[array<int, 4>{3,1, 1,-1}] = -0.20d;
	V[array<int, 4>{3,1, 3,-1}] = -0.89d;
	V[array<int, 4>{3,1,-2, 0}] = -0.97d;
	V[array<int, 4>{3,1, 0, 0}] =  2.38d;
	V[array<int, 4>{3,1,-3, 1}] =  0.89d;
	V[array<int, 4>{3,1, 1, 1}] =  1.97d;

	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
//	cout << OKQ_table[make_pair(3, 1)].format(CleanFmt) << endl;
//	cout << OKQ_table[make_pair(3, -1)].format(CleanFmt) << endl;
//

	ofstream out("iei.dat");
	for (int K1=1; K1<5; ++K1)
		for (int K2=1; K2<5; ++K2)
			for (int Q1=-K1; Q1<=K1; ++Q1)
				for (int Q2=-K2; Q2<=K2; ++Q2){
					out << K1 << " " << K2 << " " << Q1 << " " << Q2 << " ";
					out << V[array<int, 4>{K1,K2,Q1,Q2}] << endl;
				}
	out.close();

	
	for (int K1=1; K1<5; ++K1)
		for (int K2=1; K2<5; ++K2)
			for (int Q1=-K1; Q1<=K1; ++Q1)
				for (int Q2=-K2; Q2<=K2; ++Q2)
					Ham += V[array<int, 4>{K1,K2,Q1,Q2}]*OKQ_table[make_pair(K1,Q1)] * OKQ_table[make_pair(K2,Q2)];
	
	cout << Ham.format(CleanFmt) << endl;
	
	return 0;
}
