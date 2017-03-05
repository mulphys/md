/*! \file job.cc
 * \brief The ReMoDy backend 
 *
 * Executes background job.
 */
#include <string.h>
#include <cstdlib>
#include <iostream>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <math.h>
using namespace std;
#include "def.h"
#include "io.h"
#include "run.h"
#include "list.h"
#include "collection.h"
#include "model.h"
#include "species.h"
using namespace Species;
#include "domain.h"

void usage() {
	cout <<"Usage:\n\tFirst run:\n"<<Run::programname<<"\nSubsequent runs:\n"<<Run::programname<<" restart_file\n";
}
/*! \brief Strips the suffix from the filename
 * \param filename: file name including suffix
 * \return name: file name with suffix removed
 */
void parsename(char *name, char *filename) {
	char *period=strchr(filename, '.');
	int len=period==NULL?strlen(filename):period-filename;
	strncpy(name, filename, len);
	name[len]='\0';
}
/*! \brief The Main routine for a batch mode.
 * Executes a background job.
 * \param argc, argv: standard command line arguments
 * \return exit code 0 - success, 1 - failure.
 */
int	main(int argc, char *argv[]) {
	Run::init(argc,argv);//! parses command line arguments 
	Domain domain(Run::inputfile);//! creates a domain instance
	if(Run::option.xterm)printf("%c[2J",ESC);//! erases screen
	int niter=IO::getIter(Run::configfile);
	if(niter<=0) { 
		cout<<"Number of iterations "<<niter<<" should be > 0\nAborting\n";
		exit(0);
	}
	domain.run(IO::getIter(Run::configfile));//! runs the code
	char taskname[WORDLENGTH];
	if(Run::option.restart) parsename(taskname, Run::inputfile);//! retrieves taskname 
	else parsename(taskname, Run::programname);//! retrieves taskname 
	domain.save(taskname);//! saves data
	cout<<"RUN FINISHED\n";cout.flush();
	return 0;
}
