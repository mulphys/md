/*! \file view.cc
 * \brief The frontend of ReMoDy with the OpenGL Visualizer
 *
 * Executes a job and displays the computational domain 
 * in an OpenGL window.
 */
#include <iostream>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <math.h>
using namespace std;
#include "def.h"
//#include "io.h"
#include "run.h"
#include "list.h"
#include "collection.h"
#include "grid.h"
#include "container.h"
#include "model.h"
#include "species.h"
using namespace Species;
#include "domain.h"
#include "gui.h"

void usage() {
	cout <<"Usage:\n\tFirst run:\n"<<Run::programname<<"\nSubsequent runs:\n"<<Run::programname<<" restart_file\n";
}
//! The main routine
int	main(int argc, char *argv[]) {
	Run::init(argc,argv);//! parses command line arguments
//-	if(Run::inputfile[0]=='\0') {
//-		cerr<<"No input file specified\n";
//-		usage();
//-		exit(1);
//-	}
	Domain domain(Run::inputfile);//! creates a domain instance
	Gui::init(argc,argv,&domain);//! initializes the domain
	Gui::run();//! executes the job
	return 0;
}
