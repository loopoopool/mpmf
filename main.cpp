#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <array>
#include <complex>
#include "operator.h"

using namespace std;
using namespace std::complex_literals;

#ifndef D
#define D 5
#endif

complex<double> avg(double beta, Eigen::Matrix<complex<double>, D, D> O, Eigen::Vector<double, D> v){
	complex<double> tmp{0.0d}, tmp1{0.0d}, tmp2{0.0d};
	for(int i=0; i<D; ++i){
		tmp = exp(-beta*v(i));
		tmp1 += O(i,i)*tmp;
		tmp2 += tmp;
	}
	return tmp1/tmp2;
}

int main(){
	double sgn, vrcf{17.1d/120.0d};
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
			OKQ_table[make_pair(K, Q)] = pow(0.5d, 0.5d)*( sgn*top<D>(K, Q).get_matrix() + 
					top<D>(K,-Q).get_matrix() );
			OKQ_table[make_pair(K, -Q)] = 1.0i*pow(0.5d, 0.5d)*( top<D>(K, Q).get_matrix() - 
					sgn*top<D>(K,-Q).get_matrix() );
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

	// reading iei matrix
	ofstream out("iei.dat");
	for (int K1=1; K1<5; ++K1)
		for (int K2=1; K2<5; ++K2)
			for (int Q1=-K1; Q1<=K1; ++Q1)
				for (int Q2=-K2; Q2<=K2; ++Q2){
					out << K1 << " " << K2 << " " << Q1 << " " << Q2 << " ";
					out << V[array<int, 4>{K1,K2,Q1,Q2}] << endl;
				}
	out.close();


	// adding the remnant CF
	ifstream file("top_to_stevens.dat");	
	string line;
	map<pair<int, int>, double> tts_conv;
	while( getline(file, line) ){
		istringstream lineStream(line);
		vector<double> row;
		double x;
		while( lineStream >> x )
			row.push_back(x);
		tts_conv[make_pair(int(row[0]), int(row[1]))] = row[2];
	}
	file.close();
	

	srand((unsigned int) time(0));
	int nn = 6;

	Eigen::Matrix<complex<double>, D, D> tmp;
	Eigen::Matrix<complex<double>, D, D> id = Eigen::Matrix<complex<double>, D, D>::Identity();
	Eigen::Matrix<complex<double>, D, D> *occ = new Eigen::Matrix<complex<double>, D, D>[nn];

	for (int i=0; i<nn; ++i)
		occ[i] = 5.0d*Eigen::Matrix<complex<double>, D, D>::Random();
	
	map<tuple<int, int, int>, complex<double>> OKQ_mf;
	map<pair<int, int>, complex<double>> heff;

	for (int j=0; j<nn; ++j){
		tmp = occ[j] * occ[j].adjoint();
		occ[j] = tmp / tmp.trace();
	}
	
	for (int j=0; j<nn; ++j)
		for (int K=1; K<5; ++K)
			for (int Q=-K; Q<=K; ++Q){
				tmp = occ[j] * OKQ_table[make_pair(K,Q)];
				OKQ_mf[make_tuple(j,K,Q)]  = tmp.trace();
			}
	
	for(int j=0; j<nn; ++j){
		cout << "--------------------------------------------------" << endl;
		cout << "--------------------------------------------------" << endl;
		cout << " SITE : " << j << endl;
		for (int K=1; K<5; ++K){
			cout << "**************************************************" << endl;
			cout << " RANK : " << K << endl;
			for (int Q=-K; Q<=K; ++Q)
				cout << fixed << setprecision(4) << OKQ_mf[make_tuple(j,K,Q)] << "  ";
			cout << endl;
		}
	}

	int nsteps{100};
	double e0{1e10}, e{0.0d};
	double beta{2.0d}, tol{1e-6};

	// SCF loop

	for(int n=0; n<nsteps; ++n){
	
		for (int K=1; K<5; ++K)
			for (int Q=-K; Q<=K; ++Q)
				heff[make_pair(K, Q)] = 0.0d;

		for (int j=0; j<nn; ++j)
			for (int K1=1; K1<5; ++K1)
				for (int Q1=-K1; Q1<=K1; ++Q1)
					for (int K2=1; K2<5; ++K2)
						for (int Q2=-K2; Q2<=K2; ++Q2)
							heff[make_pair(K1, Q1)] += V[array<int, 4>{K1,K2,Q1,Q2}] * OKQ_mf[make_tuple(j,K2,Q2)];

		tmp = Eigen::Matrix<complex<double>, D, D>::Zero();
		e = 0.0d;
		for (int j=0; j<nn; ++j){

			for (int K=1; K<5; ++K)
				for (int Q=-K; Q<=K; ++Q)
					tmp += heff[make_pair(K, Q)] * (OKQ_mf[make_tuple(j,K,Q)]*id - 
							OKQ_table[make_pair(K,Q)]); 

			tmp += -vrcf*(OKQ_table[make_pair(4,0)]/tts_conv[make_pair(4,0)] + 
					5.0d/tts_conv[make_pair(4,4)]*OKQ_table[make_pair(4,4)]);
			
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, D, D>> eigensolver(tmp);
			if (eigensolver.info() != Eigen::Success) abort();
			e += eigensolver.eigenvalues()[0];
			tmp = Eigen::Matrix<complex<double>, D, D>::Zero();
			for (int K=1; K<5; ++K)
				for (int Q=-K; Q<=K; ++Q){
					tmp = eigensolver.eigenvectors() * OKQ_table[make_pair(K,Q)] * eigensolver.eigenvectors().adjoint();
					OKQ_mf[make_tuple(j,K,Q)] = avg(beta, tmp, eigensolver.eigenvalues());
				}
		}

		/*
		cout << "delta e : " << e - e0 << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(0,3,0)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(0,3,1)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(0,3,2)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(0,3,3)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(1,3,0)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(1,3,1)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(1,3,2)] << endl;
		cout << fixed << setprecision(4) << OKQ_mf[make_tuple(1,3,3)] << endl;
		*/
		
		if (abs(e-e0) < tol) break;
		e0 = e;
	}
	//cout << fixed << setprecision(4) << tmp.format(CleanFmt) << endl;
	//
	
	for(int j=0; j<nn; ++j){
		cout << "--------------------------------------------------" << endl;
		cout << "--------------------------------------------------" << endl;
		cout << " SITE : " << j << endl;
		for (int K=1; K<5; ++K){
			cout << "**************************************************" << endl;
			cout << " RANK : " << K << endl;
			for (int Q=-K; Q<=K; ++Q)
				cout << fixed << setprecision(4) << OKQ_mf[make_tuple(j,K,Q)] << "  ";
			cout << endl;
		}
	}
	

	
	return 0;
}
