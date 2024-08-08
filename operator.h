#ifndef OPERATOR_H
#define OPERATOR_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>


inline double fac( int &n ){
	if ( n<0 ) return 0.0d;
	else return std::tgamma( n + 1 );
}

inline double fac( int &&n ){
	return fac( n );
}



double w3j(double j1, double j2, double j3, double m1, double m2, double m3){
	if ( abs(m1+m2+m3) > 1e-9 ) return 0.0d;
	else{
		double res = fmod(j2-j1-m3, 2.0d) == 0 ? 1.0d : -1.0d;
		double tmp{1.0d}, tmp_sum{0.0d};
		int K = std::max({0.0d, j2-j3-m1, j1-j3+m2});
		int N = std::min({j1+j2-j3, j1-m1, j2+m2});
		//std::cout << " K : " << K << std::endl;
		//std::cout << " N : " << N << std::endl;
		res *= sqrt( fac(j1+j2-j3) * fac(j1-j2+j3) * fac(-j1+j2+j3) / fac(j1+j2+j3+1) );
		res *= sqrt( fac(j1-m1) * fac(j1+m1) );
		res *= sqrt( fac(j2-m2) * fac(j2+m2) );
		res *= sqrt( fac(j3-m3) * fac(j3+m3) );

		for( int k=K; k<=N; ++k ){
			tmp = k%2 == 0 ? 1.0d : -1.0d;
			//std::cout << " tmp : " << tmp << std::endl;
			tmp /= fac( k );
			tmp /= fac( j1+j2-j3-k );
			tmp /= fac( j1-m1-k );
			tmp /= fac( j2+m2-k);
			tmp /= fac( j3-j2+m1+k );
			tmp /= fac( j3-j1-m2+k );
			tmp_sum += tmp;
		}
		
		//std::cout << " tmp_sum : " << tmp_sum << std::endl;

		return res*tmp_sum;
	}
}


template <int D> class top{
	public:
		top() = default;
		top(int, int);
		top(const top &);
		void print() const;
		Eigen::Matrix<double, D, D> get_matrix() const;
	private:
		double J;
		int K, Q;
		Eigen::Matrix<double, D, D>  tkq;
};


template <int D> top<D>::top( int kk, int qq):K(kk),Q(qq){
	double sgn;
	double M1;
	J = 0.5d*(D - 1);
	for (int m1=0; m1<D; ++m1){
		M1 = m1 - J;
		sgn = fmod(J-M1, 2.0d) == 0 ? 1.0d : -1.0d;
		for (int m2=0; m2<D; ++m2)
			tkq(m1, m2) = sgn*sqrt(2.0d*K + 1.0d)*w3j(J,J,K,m2-J,-M1,Q);
	}
	
}

template <int D> Eigen::Matrix<double, D, D> top<D>::get_matrix() const{
	return tkq;
}

template <int D> void top<D>::print() const{
	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::cout << tkq.format(CleanFmt) << std::endl;
}


#endif
