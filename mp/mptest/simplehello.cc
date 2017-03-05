#include <omp.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
 {
   #pragma omp parallel  
   cout << "Hello, world.\n";
   return 0;
 }

