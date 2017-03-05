#include <stdio.h>
#include <math.h>
/*
 * Usage: hist < input > output
 * Input is a sorted list of numbers, like
 * -102.0
 *  -20.0
 *  0.1
 *  10.1
 *  ...
 * Output is the count of occurrences in each slot
 * 4
 * 2
 * 10
 * 14
 * ...
 */
#define LINELENGTH	510
const int N=100; //number of points on the histogram
const double
	xmin=-500.0, // the minimum value 
	slotwidth=10.0;

int H[N]; //the histogram

main() {
	double x=xmin;//slot coordinate
	char buf[LINELENGTH],*s;
	for (int i=0;i<N;i++)H[i]=0;
	while((s=fgets(buf,LINELENGTH,stdin))!=NULL) {
		double x;
		sscanf(s,"%lg",&x);
//		printf("z=%g ->",z);
		int i=(int)floor((x-xmin)/slotwidth);
		if(i>=0&&i<N)H[i]++;
//		printf("s=%s: x=%g, i=%d\n",s,x,i);
	}
	for(int i=0;i<N;i++) 
		printf("%d\n",H[i]);
}
