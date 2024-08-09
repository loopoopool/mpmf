#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <Eigen/Dense>

using namespace std;

int main(){
	map<pair<int, int>, double> conv_table;
	
	// source: https://doi-org/10.1016/0022-2364(90)90088-Q

	conv_table[make_pair(0, 0)] = 1.0d;
	// l = 1
	conv_table[make_pair(1, 0)] = 1.0d;
	conv_table[make_pair(1, 1)] = 1.0d;
	// l = 2
	conv_table[make_pair(2, 0)] = sqrt(6.0d);
	conv_table[make_pair(2, 1)] = sqrt(0.5d);
	conv_table[make_pair(2, 2)] = sqrt(2.0d);
	// l = 3
	conv_table[make_pair(3, 0)] = sqrt(10.0d);
	conv_table[make_pair(3, 1)] = 2.0d*sqrt(5.0d/3.0d);
	conv_table[make_pair(3, 2)] = sqrt(2.0d/3.0d);
	conv_table[make_pair(3, 3)] = 2.0d;
	// l = 4
	conv_table[make_pair(4, 0)] = 2.0*sqrt(70.0d);
	conv_table[make_pair(4, 1)] = sqrt(7.0d);
	conv_table[make_pair(4, 2)] = sqrt(14.0d);
	conv_table[make_pair(4, 3)] = 1.0d;
	conv_table[make_pair(4, 4)] = sqrt(8.0d);
	// l = 5
	conv_table[make_pair(5, 0)] = 6.0d*sqrt(14.0d);
	conv_table[make_pair(5, 1)] = 2.0d*sqrt(42.0d/5.0d);
	conv_table[make_pair(5, 2)] = sqrt(6.0d/5.0d);
	conv_table[make_pair(5, 3)] = 12.0d/sqrt(5.0d);
	conv_table[make_pair(5, 4)] = sqrt(8.0d/5.0d);
	conv_table[make_pair(5, 5)] = 4.0d;
	// l = 6
	conv_table[make_pair(6, 0)] = 4.0d*sqrt(231.0d);
	conv_table[make_pair(6, 1)] = 2.0d*sqrt(11.0d);
	conv_table[make_pair(6, 2)] = 4.0d*sqrt(22.0d/5.0d);
	conv_table[make_pair(6, 3)] = 2.0d*sqrt(22.0d/5.0d);
	conv_table[make_pair(6, 4)] = 4.0d*sqrt(11.0d/3.0d);
	conv_table[make_pair(6, 5)] = 2.0d*sqrt(2.0d/3.0d);
	conv_table[make_pair(6, 6)] = 4.0d*sqrt(2.0d);

	for (int k=1; k<=6; ++k)
		for (int q=1; q<=k; ++q)
			conv_table[make_pair(k, -q)] = conv_table[make_pair(k, q)];

	ofstream out("top_to_stevens.dat");
	for (int k=0; k<=6; ++k)
		for (int q=-k; q<=k; ++q)
			out << k << " " << q << " " << setprecision(11) << conv_table[make_pair(k,q)] << endl;	
	out.close();

	return 0;
}
