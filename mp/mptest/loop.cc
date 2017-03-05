#include <omp.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
	const int N = 1000;
	int i, a[N];
	#pragma omp parallel for
	for (i=0;i<N;i++) {
		a[i]= 2*i;
		cout << "a[" << i << "]=" << a[i] << endl;
	}
	return 0;
}
