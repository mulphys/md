#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <libxml/xmlmemory.h>
//#include <libxml/parser.h>
#include <GL/glut.h>
#include <GL/glx.h>
#include <GL/glu.h>
#include <X11/keysym.h>
#include "def.h"
#include "run.h"
#include "io.h"
using namespace std;
#include "list.h"
#include "collection.h"
#include "model.h"
#include "species.h"
using namespace Species;
#include "domain.h"
#include "gui.h"

namespace	Gui {
	Domain *domain;
	ElementDisp disp[maxelements];
	char	*configfile=(char*)"gui.cfg",
		windowname[MAXLINLEN];	
	int
		showpar[maxshowpars],
		finished,
		animation,
		nwdump,iwdump,
		attributeList[] = 
		{	GLX_RGBA, 
			GLX_RED_SIZE, 1, 
			GLX_GREEN_SIZE, 1,
			GLX_BLUE_SIZE, 1, 
			GLX_DOUBLEBUFFER, 
			GLX_DEPTH_SIZE, 1, 
			None 
		};
	REAL
		step,
		vecval[maxlineprm],
		axes[maxaxesprm],
		rgbcolor[maxcolor][3],
		wdtime,
		xo[3], //origin of coordinate system
		dx,dy,dz,
		lastx, lasty,	lastz;
	
	// Added		
	int mouseButtons[3];
	REAL zoom;
	REAL rotx;
	REAL roty;
	REAL tx;
	REAL ty;	
		
	REAL
		xmin[DIM],xmax[DIM],
		lx,ly,lz,lmin,lmax;

	enum Movement	movement;

	struct	ColorScale	colorscale;
	struct	Scene	scene;
	struct	WindowGeom	window;

	struct	BFaceList	*bface_root;

	Display *dpy;
	Window win;

	int firstR = 1;
	int *vars, *coms;

int	query_extension(char* extName) 
{
	char *p = (char *) glGetString(GL_EXTENSIONS);
	char *end = p + strlen(p);
	while (p < end) 
	{
		int n = strcspn(p, " ");
		if ((strlen(extName) == n) && (strncmp(extName, p, n) == 0))
			return GL_TRUE;
		p += (n + 1);
	}
	return GL_FALSE;
}
void	init(int argc, char* argv[], Domain *newdomain) {	

	domain=newdomain;
	
	scene.color.scheme=colorByMass;
	glutInit(&argc, argv);
	mouseButtons[0] = 0;
	mouseButtons[1] = 0;
	mouseButtons[2] = 0;
	zoom = 15.0f;
	rotx = 0.0f;
	roty = 0.001f;
	tx = 0.0f;
	ty = 0.0f;
	lastx = 0.0;
	lasty = 0.0;
	showpar[dumpWindow] = 0;
	showpar[showRun] = 1;
	showpar[showAxes] = 1;
	showpar[showSpheres] = 0;
	showpar[showNodes] = 1;
	showpar[showVariables] = 0;
	showpar[showBoundaryVertexes] = 0;
	showpar[showBoundaryVectors] = 0;
	showpar[showToolVertexes] = 0;
	showpar[showFrame] = 0;
	showpar[showBoundaryFaceCenters] = 0;
	showpar[showBoundaryFaces] = 0;
	showpar[showBoundaryGrid] = 0;
	showpar[showToolGrid] = 0;
	showpar[showCellCenters] = 0;
	showpar[showFaceCenters] = 0;
	showpar[showBonds] = 0;
	showpar[showGrid] = 0;
	finished=0;
	animation=0;
	lastz=1.0;
//-	strncpy(configfile,"gui.cfg",MAXLINLEN);//DEBUG: should be one cfg file for all
	wdtime=-1.0;
	XVisualInfo *vi;
	XSetWindowAttributes swa;
	GLXContext cx;
	sprintf(windowname,"%s",argv[0]);
	initdisp();
	dpy = XOpenDisplay(0);
	if (!dpy) ERROR("CAN'T OPEN DISPLAY");
	vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeList);
	if (!vi) ERROR("NO SUITABLE VISUAL");
	cx = glXCreateContext(dpy, vi, 0, GL_TRUE);
	swa.colormap = XCreateColormap
	(	dpy, 
		RootWindow(dpy, vi->screen),
		vi->visual, AllocNone
	);
	swa.border_pixel = 0;
	swa.event_mask = ExposureMask | StructureNotifyMask | KeyPressMask |
	ButtonPressMask | ButtonMotionMask;
	if(Run::option.debug)
		printf
		(	"Open window: width=%d, height=%d\n",
			window.width,window.height
		);
	if(Run::option.debug)
		printf("Window: width=%d, height=%d\n",window.width, window.height);
	win=XCreateWindow
	(	dpy, 
		RootWindow(dpy, vi->screen), 
		0, 0, window.width, window.height,
		0, vi->depth, InputOutput, vi->visual,
		CWBorderPixel|CWColormap|CWEventMask, &swa
	);
	XStoreName(dpy, win, windowname);
	XMapWindow(dpy, win);
	glXMakeCurrent(dpy, win, cx);
	/* check for the polygon offset extension */
#ifndef GL_VERSION_1_1
	if (!query_extension("GL_EXT_polygon_offset"))
		error(argv[0], "polygon_offset extension is not available");
#else
	(void) query_extension;
#endif
	Materials(argc, argv);
	helpDisplay();
}
void	helpDisplay() {
	printf("MOUSE CONTROL:\n");
	printf("left button & drag - rotate\n");
	printf("right button & drag - translate\n");
	printf("KEYSTROKE COMMANDS:\n");
	printf("?     - show help\n");
	printf("a     - show/hide axes\n");
	printf("c     - switch color scheme\n");
	printf("f     - show/hide frame\n");
	printf("G(g)  - show/hide (boundary) grid\n");
	printf("i     - initialize domains\n");
	printf("l     - load configuration from %s\n",configfile);
	printf("n     - show/hide nodes\n");
	printf("o     - toggle output\n");
	printf("q     - quit without saving\n");
	printf("r     - run\n");
	printf("s     - advance one iteration step\n");
	printf("t(T)  - toggle pores nodes\n");
	printf("w     - write data at this time step into %s-<number>.dat.gz\n",Run::outputname);
	printf("+     - move forward along the arrow-axis\n");
	printf("-     - move backward along the arrow-axis\n");
	printf(":     - command mode\n esc: exit\n");
}
int	readparam(//TO BE DEPRICATED
	char	*s,
	char	*param[],
	int	maxparam,
	char	*filename,
	int	val[]
)
{	int	i;
	for (i=0; i<maxparam; i++)
	{	int lparam;
		if(param[i]==NULL) continue;
		lparam=strlen(param[i]);
		if 
		(	lparam>0 &&
			(*(s+lparam)==' ' || *(s+lparam)=='=')
			&& memcmp(s,param[i],lparam)==0
		)	break;
	}
	if (i<maxparam)
	{	if ((s=strchr(s,'='))==NULL)
		{	fprintf
			(	stderr,
				"WARNING: wrong format at '%s' in file %s\n\t-> IGNORED\n",
				s,filename
			);
			return -1;
		}
		while(isspace(*++s));
		sscanf(s,"%d",&val[i]);
	}
	return i;
}
int	readparam(//TO BE DEPRICATED
	char	*s,
	char	*param[],
	int	maxparam,
	char	*filename,
	REAL	val[]
)
{	int	i;
	for (i=0; i<maxparam; i++)
	{	int lparam=strlen(param[i]);
		if 
		(	(*(s+lparam)==' ' || *(s+lparam)=='=')
			&& memcmp(s,param[i],lparam)==0
		)	break;
	}
	if (i<maxparam)
	{	if ((s=strchr(s,'='))==NULL)
		{	fprintf
			(	stderr,
				"WARNING: wrong format at '%s' in file %s\n\t-> IGNORED\n",
				s,filename
			);
			return -1;
		}
		while(isspace(*++s));
		sscanf(s,"%g",&val[i]);
	}
	return i;
}
void	readconf() {	
	using namespace IO;
	int	i;
	char
		*colorname[maxcolor],
		*showparnam[maxshowpars],
		*pointparam[maxpointprm],
		*lineparam[maxlineprm],
		*axesprmnam[maxaxesprm],
//-		*poreprmnam[maxporeprm],
		*filename=configfile,
		buf[MAXLINLEN],*s;
	REAL	a;
	FILE	*file;
	colorname [red      ]=strdup("color.red");
	colorname [green    ]=strdup("color.green");
	colorname [blue     ]=strdup("color.blue");
	colorname [skyblue  ]=strdup("color.skyblue");
	colorname [brown    ]=strdup("color.brown");
	colorname [magenta  ]=strdup("color.magenta");
	for (int i=magenta+1;i<maxcolor;i++)
	{	char	name[MAXLINLEN];
		sprintf(name,"color%d",i);
		colorname[i]=strdup(name);
	}
//BEGIN DEPRECATION: OBSOLETE: used for gui.cfg input. Will be deprecated.
	for (int i=0;i<maxshowpars;i++)
		showparnam[i]=NULL;
	showparnam[showSpheres]=strdup("spheres"    );
	pointparam[variable ]=strdup("variable"     );
	pointparam[size     ]=strdup("size"         );
	pointparam[sizevar  ]=strdup("sizevar"      );
	pointparam[massvar  ]=strdup("massvar"      );
	pointparam[vmin     ]=strdup("value.min"    );
	pointparam[vmax     ]=strdup("value.max"    );
	lineparam [thickness]=strdup("thickness"    );
	lineparam [length   ]=strdup("length"       );
	axesprmnam[axeswidth]=strdup("width"        );
	axesprmnam[xaxislength]=strdup("length.x"   );
	axesprmnam[yaxislength]=strdup("length.y"   );
	axesprmnam[zaxislength]=strdup("length.z"   );
	axesprmnam[arrowidth  ]=strdup("arrow.width");
	axesprmnam[arrowheight]=strdup("arrow.length");
	for (i=0; i<maxcolor; i++)
	{	pointparam[i]=strdup(colorname[i]);
		lineparam [i]=strdup(colorname[i]);
//-		poreprmnam [i]=strdup(colorname[i]);
	}
//-	if ((file=fopen(filename,"r"))==NULL)
//-	{	fprintf
//-		(	stderr,
//-			"CAN'T OPEN CONFIGURATION FILE %s: using default settings\n",
//-			filename
//-		);return;
//-	}
	step=1.0;//default translation step
	int colorscheme=0;//color scheme
	scene.mesh.node[size]=1.0;
	scene.mesh.node[vmin]=0.0;
	scene.mesh.node[vmax]=1.0;
	scene.mesh.node[red]=0.0;
	scene.mesh.node[green]=1.0;
	scene.mesh.node[blue]=0.0;
	scene.mesh.line[thickness]=1.0;
	scene.mesh.line[red]=1.0;
	scene.mesh.line[green]=0.0;
	scene.mesh.line[blue]=0.0;
	scene.frame.line[thickness]=1.0;
	scene.frame.line[red]=1.0;
	scene.frame.line[green]=1.0;
	scene.frame.line[blue]=0.0;
//-	if(Run::option.verbose)
//-		printf("Reading configuration from %s ...",filename);fflush(stdout);
//-	while (fgets(buf,MAXLINLEN,file)!=NULL) {
//-		for (s=buf;isspace(*s);s++);
//-		if (memcmp(s,"Translation.step",16)==0) {
//-			if ((s=strchr(s,'='))==NULL) {
//-				fprintf(stderr,"MISSING '=' in Step difinition in %s\n",filename);
//-				continue;
//-			}
//-			while(isspace(*++s));
//-			sscanf(s,"%g",&step);
//-			dx=step*lx;
//-			dy=step*ly;
//-			dz=step*lz;
//-		}
//-		else
//-		if (memcmp(s,"Display.",8)==0) {
//-			if (readparam(s+8,showparnam,maxshowpars,filename,showpar)==-1)
//-				continue;
//-		}
//-		else
//-		if (memcmp(s,"Vector.",7)==0) {
//-			if (readparam(s+7,lineparam,maxlineprm,filename,vecval)==-1)
//-				continue;
//-		}
//-		else
//-		if (memcmp(s,"Grid.",5)==0) {
//-			s+=5;
//-			if (memcmp(s,"node.",5)==0)
//-			{	if 
//-				(	(i=readparam(s+5,pointparam,maxpointprm,filename,scene.node))==-1
//-				)	continue;
//-			}
//-			else
//-			if (memcmp(s,"line.",5)==0) {
//-				if 
//-				(	(i=readparam(s+5,lineparam,maxlineprm,filename,scene.line))==-1
//-				)	continue;
//-			}
//-		}
//-		else
//-		if (memcmp(s,"Axes.",5)==0)
//-		{	if 
//-			(	(i=readparam(s+5,axesprmnam,maxaxesprm,filename,axes))==-1
//-			)	continue;
//-		}
//-	}
//-	fclose(file);
//-	if(Run::option.verbose)printf(" done\n");fflush(stdout);
//END DEPRECATION
	for(int i=0;i<DIM;i++) {
		xmin[i]=domain->minBound(i);
		xmax[i]=domain->maxBound(i);
	}
	//Reading XML input from Run::configfile
	char *doctype=(char*)DOCTYPE;
	char *inpfilename=Run::configfile;
//-	IO::getTime(inpfilename);
	//READ DATA INPUT FILE
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	if(Run::option.verbose)printf("Domain::inpfilename=%s\n",inpfilename);
	xmlDocPtr doc;
	xmlNodePtr root;
	//Corrdinates of basic element types
	if (Run::option.verbose|Run::option.debug)
	{	printf
		(	"Reading geometry from %s\n",
			inpfilename
		);fflush(stdout);
	}
	doc = xmlParseFile(inpfilename);
	if (doc == NULL ) 
	{	fprintf(stderr,"XML parser failed in %s\n",inpfilename);
		exit(1);
	}
	root = xmlDocGetRootElement(doc);
	if (root == NULL) 
	{	fprintf(stderr,"Empty document: %s\n",inpfilename);
		xmlFreeDoc(doc);
		return;
	}
	if (xmlStrcmp(root->name, (const xmlChar *) doctype)) 
	{	fprintf
		(	stderr,
			"document of the wrong type, root node != %s in %s\n",
			doctype,inpfilename
		);
		xmlFreeDoc(doc);
		return;
	}
	for 
	(	xmlNodePtr cur = root->xmlChildrenNode; 
		cur != NULL; cur = cur->next
	) 
	{	//SCAN DESCRIPTION AND FIND THIS DOMAIN
//-		if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))) {
//-			char name[MAXLINLEN];
//-			if(getCharAttr(cur,"name",name)==0) {
//-				fprintf
//-				(	stderr,"CAN'T GET %s FOR %s IN %s\nAborting\n",
//-					"name","domain",inpfilename
//-				); exit(1);
//-			}
//-			if(Run::option.verbose) printf("Domain: %s\n",name);
//-			if(parseWord(doc,cur,"type",word)==0) {
//-				fprintf
//-				(	stderr,"Tag \"type\" is missing in %s\nAborting\n",
//-					inpfilename
//-				); exit(1);
//-			}
//-			if(Run::option.verbose){printf("Domain type: %s\n",word);fflush(stdout);}
//-			for 
//-			(	xmlNodePtr next=cur->xmlChildrenNode;
//-				next != NULL;
//-				next=next->next
//-			) {
//-				if ((!xmlStrcmp(next->name, (const xmlChar *)"geometry"))) {
//-					if(parseWord(doc,next,"bounds",word)!=0) {
//-						sscanf
//-						(	word,"%g %g %g %g %g %g",
//-							xmin,xmax,xmin+1,xmax+1,xmin+2,xmax+2
//-						);
//-						if(Run::option.verbose)
//-						printf
//-						(	"\tbounds: (%g:%g), (%g:%g), (%g:%g)\n",
//-							xmin[0],xmax[0],xmin[1],xmax[1],xmin[2],xmax[2]
//-						);
//-						//Check domain bounds
//-						bool adjusted=false;
//-						for(int i=0;i<DIM;i++) {
//-							if(xmin[i]>domain->minBound(i)) {
//-								xmin[i]=domain->minBound(i);
//-								adjusted=true;
//-							} else domain->setMinBound(i,xmin[i]);
//-							if(xmax[i]<domain->maxBound(i)) {
//-								xmax[i]=domain->maxBound(i);
//-								adjusted=true;
//-							} else domain->setMaxBound(i,xmax[i]);
//-						}
//-						if(adjusted) 
//-							printf
//-							(	"\tDomain bounds adjusted: (%g:%g), (%g:%g), (%g:%g)\n",
//-								xmin[0],xmax[0],xmin[1],xmax[1],xmin[2],xmax[2]
//-							);
//-					}
//-					if(parseWord(doc,next,"boundary",word)!=0) {
//-						if(strcmp(word,"periodic")==0) domain->BoundaryType(periodicBoundary);
//-						else
//-						if(strcmp(word,"elastic")==0) domain->BoundaryType(elasticBoundary);
//-						else domain->BoundaryType(periodicBoundary);
//-					}
//-				}//END geometry
//-			}//END next
//-		}//END if domain
		//READ GUI SECTION
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"gui"))) {//GUI PARAMETERS:
			for(xmlNodePtr gui = cur->xmlChildrenNode; gui != NULL; gui = gui->next) {
				// Find the translation
				if((!xmlStrcmp(gui->name, (const xmlChar *)"translation"))) {
					if( parseWord(doc, gui, (char*)"step", word) != 0 ) {
						REAL step = (REAL)atof(word);
						if(Run::option.verbose) cout << "Translation step = " << step << endl;
					}
				}
				// Find the vector
				if((!xmlStrcmp(gui->name, (const xmlChar *)"vector"))) {	
					if(Run::option.verbose) cout<<"Vector parameters:\n";
					if( parseWord(doc, gui, (char*)"length", word) != 0 ) {
						REAL veclength = (REAL)atof(word);
						if(Run::option.verbose) cout << "\tlength = " << veclength << endl;
						vecval[length]=veclength;
					}
					
					if( parseWord(doc, gui, (char*)"thickness", word) != 0 ) {
						REAL vecthickness = (REAL)atof(word);
						if(Run::option.verbose) cout << "\tthickness = " << vecthickness << endl;
						vecval[thickness]=vecthickness;
					}
					for(xmlNodePtr color = gui->xmlChildrenNode; color != NULL; color = color->next) {
						if ((!xmlStrcmp(color->name, (const xmlChar *)"color"))) {
							if( parseWord(doc, color, (char*)"red", word) != 0 ) {
								REAL vectorRed = (REAL)atof(word);
								if(Run::option.verbose) 
									cout << "\tColor Red = " << vectorRed << endl;
								vecval[red]=vectorRed;
							}
					
							if( parseWord(doc, color, (char*)"green", word) != 0 ) {
								REAL vectorGreen = (REAL)atof(word);
								if(Run::option.verbose) 
									cout << "\tColor Green = " << vectorGreen << endl;
								vecval[green]=vectorGreen;
							}
					
							if( parseWord(doc, color, (char*)"blue", word) != 0 ) 
							{
								REAL vectorBlue = (REAL)atof(word);
								if(Run::option.verbose) 
									cout << "\tColor Blue = " << vectorBlue << endl;
								vecval[blue]=vectorBlue;
							}							
						}// end color
					}// End color loop
				}// end vector
				if((!xmlStrcmp(gui->name, (const xmlChar *)"frame"))) {
					if(Run::option.verbose)cout<<"Frame parameters\n";
					for( xmlNodePtr fr = gui->xmlChildrenNode; 
						fr != NULL; fr = fr->next 
					) {
						if ((!xmlStrcmp(fr->name, (const xmlChar *)"line"))) {
							if(Run::option.verbose)cout<<"\tLine\n";
							if( parseWord(doc, fr, (char*)"thickness", word) != 0 ) {
								REAL lineThickness = (REAL)atof(word);
								if(Run::option.verbose) cout << "\t\tThickness = " << lineThickness << endl;
								scene.frame.line[thickness]=lineThickness;
							}
							for(xmlNodePtr color = fr->xmlChildrenNode; 
								color != NULL; color = color->next) {
								if ((!xmlStrcmp(color->name, (const xmlChar *)"color"))) {
									if(Run::option.verbose)cout<<"\t\tColor\n";
									if( parseWord(doc, color, (char*)"red", word) != 0 ) {
										REAL lineRed = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tRed = " << lineRed << endl;
										scene.frame.line[red]=lineRed;
									}
									if( parseWord(doc, color, (char*)"green", word) != 0 ) {
										REAL lineGreen = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tGreen = " << lineGreen << endl;
										scene.frame.line[green]=lineGreen;
									}
									if( parseWord(doc, color, (char*)"blue", word) != 0 ) {
										REAL lineBlue = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tBlue = " << lineBlue << endl;
										scene.frame.line[blue]=lineBlue;
									}
									REAL rgbcolor[3];
									rgbcolor[0]=scene.frame.line[red];
									rgbcolor[1]=scene.frame.line[green];
									rgbcolor[2]=scene.frame.line[blue];
									setElementColor(frame,rgbcolor);
								}// endif(color)
							}// end for(color)
						}// end if(line)
					}//end for
				}//end if(frame)
				if((!xmlStrcmp(gui->name, (const xmlChar *)"mesh"))) {
					if(Run::option.verbose)cout<<"Mesh parameters\n";
					for( xmlNodePtr fr = gui->xmlChildrenNode; 
						fr != NULL; fr = fr->next 
					) {
						if ((!xmlStrcmp(fr->name, (const xmlChar *)"node"))) {
							if(Run::option.verbose)cout<<"\tNode\n";
							if( parseWord(doc, fr, (char*)"type", word) != 0 ) {
								if(Run::option.verbose) cout << "\t\tType = " << word << endl;
								showpar[showSpheres] = 0;
								if(strncmp(word,"sphere",6)==0) showpar[showSpheres]=1;
							}
							for(xmlNodePtr colors = fr->xmlChildrenNode; 
								colors != NULL; colors = colors->next) {
								if ((!xmlStrcmp(colors->name, (const xmlChar *)"colorscheme"))) {
									if(Run::option.verbose)cout<<"\t\tColor Scheme\n";
									if( parseWord(doc, colors, (char*)"variable", word) != 0 ) {
										if(Run::option.verbose) cout << "\t\tVariable = " << word << endl;
										colorscheme=0;
										if(strcmp(word,"type")==0) colorscheme=(int)colorByType;
										else if(strcmp(word,"mass")==0) colorscheme=(int)colorByMass;
									}
									scene.color.minvalue=0.0;
									if( parseWord(doc, colors, (char*)"minvalue", word) != 0 ) {
										if(Run::option.verbose) cout << "\t\tMin Value = " << word << endl;
										REAL minvalue=(REAL)atof(word);
										scene.color.minvalue=minvalue;
									}
									scene.color.maxvalue=1.0;
									if( parseWord(doc, colors, (char*)"maxvalue", word) != 0 ) {
										if(Run::option.verbose) cout << "\t\tMax Value = " << word << endl;
										REAL maxvalue=(REAL)atof(word);
										scene.color.maxvalue=maxvalue;
									}
								}
							}
							if( parseWord(doc, fr, (char*)"size", word) != 0 ) {
								REAL nodeSize = (REAL)atof(word);
								if(Run::option.verbose) cout << "\t\tSize = " << nodeSize << endl;
								scene.mesh.node[size]=nodeSize;
							}
							for(xmlNodePtr color = fr->xmlChildrenNode; color != NULL; color = color->next) {
								if ((!xmlStrcmp(color->name, (const xmlChar *)"color"))) {
									if(Run::option.verbose)cout<<"\t\tColor\n";
									if( parseWord(doc, color, (char*)"red", word) != 0 ) {
										REAL lineRed = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tRed = " << lineRed << endl;
										scene.mesh.node[red]=lineRed;
									}
									if( parseWord(doc, color, (char*)"green", word) != 0 ) {
										REAL lineGreen = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tGreen = " << lineGreen << endl;
										scene.mesh.node[green]=lineGreen;
									}
									if( parseWord(doc, color, (char*)"blue", word) != 0 ) {
										REAL lineBlue = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tBlue = " << lineBlue << endl;
										scene.mesh.node[blue]=lineBlue;
									}
									REAL rgbcolor[3];
									rgbcolor[0]=scene.mesh.node[red];
									rgbcolor[1]=scene.mesh.node[green];
									rgbcolor[2]=scene.mesh.node[blue];
									setElementColor(nodes,rgbcolor);
								}// endif(color)
							}// end for(color)
						}// end if(line)
						if ((!xmlStrcmp(fr->name, (const xmlChar *)"line"))) {
							if(Run::option.verbose)cout<<"\tLine\n";
							if( parseWord(doc, fr, (char*)"thickness", word) != 0 ) {
								REAL lineThickness = (REAL)atof(word);
								if(Run::option.verbose) cout << "\t\tThickness = " << lineThickness << endl;
								scene.mesh.line[thickness]=lineThickness;
							}
							for(xmlNodePtr color = fr->xmlChildrenNode; 
								color != NULL; color = color->next) {
								if ((!xmlStrcmp(color->name, (const xmlChar *)"color"))) {
									if(Run::option.verbose)cout<<"\t\tColor\n";
									if( parseWord(doc, color, (char*)"red", word) != 0 ) {
										REAL lineRed = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tRed = " << lineRed << endl;
										scene.mesh.line[red]=lineRed;
									}
									if( parseWord(doc, color, (char*)"green", word) != 0 ) {
										REAL lineGreen = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tGreen = " << lineGreen << endl;
										scene.mesh.line[green]=lineGreen;
									}
									if( parseWord(doc, color, (char*)"blue", word) != 0 ) {
										REAL lineBlue = (REAL)atof(word);
										if(Run::option.verbose) cout << "\t\t\tBlue = " << lineBlue << endl;
										scene.mesh.line[blue]=lineBlue;
									}
									REAL rgbcolor[3];
									rgbcolor[0]=scene.mesh.line[red];
									rgbcolor[1]=scene.mesh.line[green];
									rgbcolor[2]=scene.mesh.line[blue];
									setElementColor(mesh,rgbcolor);
								}// endif(color)
							}// end for(color)
						}// end if(line)
					}//end for
				}//end mesh
			}//end for
		}// END GUI
	}//END cur
	xmlFreeDoc(doc);
	dx=step*lx;
	dy=step*ly;
	dz=step*lz;
	scene.color.scheme=(ColorSchemes) colorscheme;
	if(Run::option.verbose){printf("File %s parsed\n",inpfilename);fflush(stdout);}
}
void	refresh() {
//-	if
//-	(	showpar[showBoundaryFaces]||showpar[showBoundaryGrid]
//-	)
//-	{
//-///		if(domain->getNoVariables(boundary_faces)==0)
//-///		{	domain->deleteRing(bface_root);
//-///			domain->createBoundaryFaceList(bface_root);
//-///		}
//-///		else bface_root=domain->bface_root;
//-	}
}
void	initdisp() {
	if(Run::option.debug){printf("Initializing display\n");fflush(stdout);}
	rgbcolor[red    ][red      ]=1.0;
	rgbcolor[red    ][green    ]=0.0;
	rgbcolor[red    ][blue     ]=0.0;
	rgbcolor[green  ][red      ]=0.0;
	rgbcolor[green  ][green    ]=1.0;
	rgbcolor[green  ][blue     ]=0.0;
	rgbcolor[blue   ][red      ]=0.0;
	rgbcolor[blue   ][green    ]=0.0;
	rgbcolor[blue   ][blue     ]=1.0;
	rgbcolor[skyblue][red      ]=0.0;
	rgbcolor[skyblue][green    ]=1.0;
	rgbcolor[skyblue][blue     ]=1.0;
	rgbcolor[brown  ][red      ]=1.0;
	rgbcolor[brown  ][green    ]=1.0;
	rgbcolor[brown  ][blue     ]=0.0;
	rgbcolor[magenta][red      ]=1.0;
	rgbcolor[magenta][green    ]=0.0;
	rgbcolor[magenta][blue     ]=1.0;
	rgbcolor[yellow][red      ]=1.0;
	rgbcolor[yellow][green    ]=1.0;
	rgbcolor[yellow][blue     ]=0.0;
	colorscale.relative=0;
	if (maxcolor>6)
	{	int
			ncolor=(int)pow((double)maxcolor,1./3.),
			icolor=5;
		for (int ired  =1; ired  <ncolor; ired++   )
		for (int igreen=1; igreen<ncolor; igreen++ )
		for (int iblue =1; iblue <ncolor; iblue++  )
		{	if (++icolor>=maxcolor) break;
			rgbcolor[icolor][ired  ]=(REAL)ired  /(REAL)ncolor;
			rgbcolor[icolor][igreen]=(REAL)igreen/(REAL)ncolor;
			rgbcolor[icolor][iblue ]=(REAL)iblue /(REAL)ncolor;			
		}
	}
	vecval[length   ]=.2;
	vecval[thickness]=1.5;
	vecval[red      ]=0.5;
	vecval[green    ]=1.0;
	vecval[blue     ]=1.0;
	axes[axeswidth  ]=5.0;
	axes[xaxislength]=0.1;
	axes[yaxislength]=0.1;
	axes[zaxislength]=0.1;
	axes[arrowheight]=0.1;
	axes[arrowidth  ]=0.03;
//	scene.ivar=-1;
	scene.frame.line[thickness]=1.0;
	scene.frame.line[red      ]=1.0;
	scene.frame.line[green    ]=0.0;
	scene.frame.line[blue     ]=0.0;
	scene.mesh.node[size ]=1.0;
	scene.mesh.node[red  ]=1.;
	scene.mesh.node[green]=1.;
	scene.mesh.node[blue ]=0.;
	scene.mesh.line[thickness]=1.0;
	scene.mesh.line[red      ]=1.0;
	scene.mesh.line[green    ]=0.0;
	scene.mesh.line[blue     ]=0.0;
	for (int eltype=0; eltype<maxelements; eltype++)
		setElementColor(eltype,rgbcolor[(eltype)%maxcolor]);
	setElementColor(nodes,rgbcolor[green]);
	setElementColor(edges,rgbcolor[red]);
	setElementColor(faces,rgbcolor[brown]);
	setElementColor(cells,rgbcolor[skyblue]);
	setElementColor(frame,rgbcolor[yellow]);
	setElementColor(mesh,rgbcolor[red]);
	readconf();
	//If domain is degenerate, assign xmin/xmax from Gui:
	if(domain->Molecules()->number()==0) {
		for(int i=0;i<DIM;i++) {
			domain->setMinBound(i,xmin[i]);
			domain->setMaxBound(i,xmax[i]);
		}
	}
	if(Run::option.debug)
	{	printf
		(	"Grid limits: xmin={%g,%g,%g}, xmax={%g,%g,%g}\n",
			xmin[0],xmin[1],xmin[2],xmax[0],xmax[1],xmax[2]
		); fflush(stdout);
	}
//	domain->setDisplay();// lighting=1;

	lmin= LARGE;
	lmax=-LARGE;
	lx=xmax[0]-xmin[0];
	ly=xmax[1]-xmin[1];
	lz=xmax[2]-xmin[2];
	if (lx>lmax) lmax=lx;
	if (ly>lmax) lmax=ly;
	if (lz>lmax) lmax=lz;
	if (lx<lmin) lmin=lx;
	if (ly<lmin) lmin=ly;
	if (lz<lmin) lmin=lz;
	xo[0]=xmin[0]+.5*lx;
	xo[1]=xmin[1]+.5*ly;
	xo[2]=xmin[2]+.5*lz;
lx=ly=lz=lmax;///DDD
	dx=step*lx;
	dy=step*ly;
	dz=step*lz;
	zoom=1.4*lz;
	if (Run::option.debug|Run::option.verbose)
		printf("Origin of coordinate axes: %g\t%g\t%g\n",xo[0],xo[1],xo[2]);
	{	REAL	window_size;
		window_size=1.5*WINDOW_SIZE/sqrt(lx*lx+ly*ly);
		window.width=(int)(window_size*lx); 
		if (window.width>MAX_WINDOW_WIDTH)window.width=MAX_WINDOW_WIDTH; 
		window.height=(int)(window_size*ly);
		if (window.height>MAX_WINDOW_HEIGHT)window.height=MAX_WINDOW_HEIGHT; 
	}
//-	domain->setBoundary();
}
/*  Initialize material property and depth buffer.
 */
void	Materials(int argc, char *argv[]) {
	GLfloat mat[4] = { 0.8, 0.7, 0.5, 1.0 },
		shine=0.6;
	if(Run::option.debug)
	{	printf("Materials: argc=%d",argc);fflush(stdout);
		for (int i=0; i<argc; i++)
		printf(" '%s'",argv[i]);
		printf("\n");fflush(stdout);
	}
//	InitMaterials();
	GLfloat ambient[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat position[] = { 0.0, 3.0, 3.0, 0.0 };
	GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat local_view[] = { 0.0 };
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);
	glFrontFace (GL_CW);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearAccum(0.0, 0.0, 0.0, 0.0);

	glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
	glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
	glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
	glTranslatef(-xo[0], -xo[1], -xo[2]-3*lz);
}
//- void	showGridElements (
//- 	int	ne,//number of elements to show
//- 	double	*X,//coordinates
//- 	REAL	color[]
//- ) {
//- 	int	i;
//- 	glPointSize(scene.node[size]);
//- 	glColor3f(color[red], color[green], color[blue]);
//- 	glDisable(GL_LIGHTING);
//- 	glBegin(GL_POINTS);
//- 	for (i=0; i<ne; i++) {
//- 		double	*x=X+DIM*i;
//- 		glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//- 	}
//- 	glEnd();
//- 	glEnable(GL_LIGHTING);
//- }
void	showVector
(
	double	*x,//coordinates
	double	*v //vector
)
{
	if (vecval[length]*(fabs(v[0])+fabs(v[1])+fabs(v[2]))>1.e6)
	{	fprintf
		(	stderr,
			"VALUE AT (%g,%g,%g) TOO LARGE (not displayed)\n",
			x[0],x[1],x[2]
		);
		return;
	}
	glBegin(GL_LINES);
	glVertex3f
	(	(GLfloat)x[0],
		(GLfloat)x[1],
		(GLfloat)x[2]
	);
	glVertex3f
	(	(GLfloat)(x[0]+vecval[length]*v[0]),
		(GLfloat)(x[1]+vecval[length]*v[1]),
		(GLfloat)(x[2]+vecval[length]*v[2])
	);
	glEnd();
}
//+ template <class Element>
//+ void	showDVectors
//+ (	int loc,
//+ 	Element *root
//+ )
//+ {
//+ 	if (root!=NULL)
//+ 	{	Element	*e=root;
//+ 		do
//+ 		{	double
//+ 				*x=e->x,
//+ 				*v=e->var+loc;
//+ 			showVector(x,v);
//+ 			e=e->next;
//+ 		}	while(e!=root);
//+ 	}
//+ }
//+ void	showVectors
//+ (
//+ 	int	ivar,
//+ 	Domain	*domain
//+ )
//+ {
//+ 	Variable	*var=domain->variable+ivar;
//+ 	int
//+ 		type=var->type,
//+ 		size=var->size;
//+ 	double
//+ 		*X=domain->coordinates[type].val,//coordinates
//+ 		*x,*v,*u;
//+ 	glLineWidth (vecval[thickness]);
//+ 	glDisable(GL_LIGHTING);
//+ 	glColor3f(vecval[red], vecval[green], vecval[blue]);
//+ 	u=var->val; //scene.V+scene.var[scene.ivec].i*scene.n*DIM;
//+ 	if (domain->type==dynamic)
//+ 	{	int	loc=var->loc;
//+ 		switch(type)
//+ 		{	case nodes:
//+ 			showDVectors(loc,domain->dnode_root);
//+ 			break;
//+ 			case cells:
//+ 	///		showDVectors(loc,domain->dcell_root);
//+ 	///		break;
//+ 			default:
//+ 			fprintf(stderr,"Can't display vector variable type %d",type);
//+ 		}
//+ 	}
//+ 	else
//+ 		for (x=X,v=u;x-X<size*DIM;x+=DIM,v+=DIM)
//+ 			showVector(x,v);
//+ 	glEnable(GL_LIGHTING);
//+ }
void	showParticle (
	double	vmn,
	double	vmx,
	double	val,
	double	*x
) {
		REAL	r,g,b,
			vav=.5*(vmn+vmx),
			alow =PI/(2.*(vav-vmn)),
			ahigh=PI/(2.*(vmx-vav)),
			a;
		if (val<vmn)
		{	r=g=0.0; b=1.0; }
		else
		if (val>vmx)
		{	r=1.0; g=b=0.0; }
		else
		if (val<vav)
		{	a=(REAL)(val-vmn)*alow;
			b=cos(a); g=sin(a);
			r=0.0;
		}
		else
		{	a=(REAL)(val-vav)*ahigh;
			g=cos(a); r=sin(a);
			b=0.0;
		}
		glColor3f(r,g,b);
		glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
}
void	showSphere (
	double	vmn,
	double	vmx,
	double	val,
	double	rad,
	double	*x
)
{
		GLfloat mat[4],/// = { 0.6, 0.8, 0.3, 1.0 },
				shine=0.6;
		REAL	r,g,b,
			vav=.5*(vmn+vmx),
			alow =PI/(2.*(vav-vmn)),
			ahigh=PI/(2.*(vmx-vav)),
			a;
		if (val<vmn)
		{	r=g=0.0; b=1.0; }
		else
		if (val>vmx)
		{	r=1.0; g=b=0.0; }
		else
		if (val<vav)
		{	a=(REAL)(val-vmn)*alow;
			b=cos(a); g=sin(a);
			if(g<0.0)g=0.0;
			if(b<0.0)b=0.0;
			r=0.0;
		}
		else
		{	a=(REAL)(val-vav)*ahigh;
			g=cos(a); r=sin(a);
			if(g<0.0)g=0.0;
			if(r<0.0)r=0.0;
			b=0.0;
		}
		
		// ADDED
		glPushMatrix();
		
		mat[0]=r;mat[1]=g;mat[2]=b;mat[3]=1.0;
		glEnable(GL_LIGHTING);
		glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
		glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
		glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
		glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
		glTranslatef (x[0], x[1], x[2]); 
		glutSolidSphere (rad, 15, 15);
		
		//glTranslatef (-x[0],-x[1],-x[2]);
		glPopMatrix();
}
//+ template <class Element>
//+ void	showDScalar
//+ (	int loc,
//+ 	Element *root,
//+ 	double	vmn,
//+ 	double	vmx
//+ )
//+ {	if (root==NULL)return;
//+ 	{	Element	*e=root;
//+ 		glBegin(GL_POINTS);
//+ 		do
//+ 		{	double value=e->var[loc];
//+ 			showPore(vmn,vmx,value,e->x);
//+ 			e=e->next;
//+ 		}	while(e!=root);
//+ 		glEnd();
//+ 	}
//+ }
//+ void	showScalar
//+ (
//+ 	int	ivar,
//+ 	int	icomp,
//+ 	double	vmn,
//+ 	double	vmx,
//+ 	Domain	*domain
//+ )
//+ {
//+ 	Variable	*var=domain->variable+ivar;
//+ 	char	*name=var->name;
//+ 	int
//+ 		type=var->type,
//+ 		dim=var->dimension,
//+ 		size=var->size,
//+ 		jcomp=icomp;
//+ 	double
//+ 		*X=domain->coordinates[type].val,
//+ 		*value=var->val;
//+ 	if (icomp>DIM)
//+ 	{	fprintf
//+ 		(	stderr,"Can't display variable component %d for variable %s\n",
//+ 			icomp,name
//+ 		);
//+ 		return;
//+ 	}
//+ 	if (dim==DIM&&icomp<0)
//+ 	{	fprintf
//+ 		(	stderr,"A vector component should first be selected for variable %s\n",
//+ 			name
//+ 		);
//+ 		return;
//+ 	}
//+ 	if (icomp>=dim)
//+ 	{	fprintf
//+ 		(	stderr,
//+ 			"Can't display component %d of variable %s of dim=%d\n",
//+ 			icomp,name,dim
//+ 		);
//+ 		return;
//+ 	}
//+ 	if (icomp<0) jcomp=0;
//+ 	if (domain->type==dynamic)
//+ 	{	int	loc=var->loc;
//+ 		switch(type)
//+ 		{	case nodes:
//+ 				showDScalar(loc,domain->dnode_root,vmn,vmx);
//+ 				break;
//+ 			case cells:
//+ 			///	showDScalar(loc,domain->dcell_root,vmn,vmx);
//+ 			///	break;
//+ 			default:
//+ 			fprintf(stderr,"Can't display scalar variable type %d",type);
//+ 		}
//+ 	}
//+ 	else
//+ 	{
//+ 		glBegin(GL_POINTS);
//+ 		for (int i=0; i<size; i++)
//+ 			showPore(vmn,vmx,value[i*dim+jcomp],X+DIM*i);
//+ 		glEnd();
//+ 	 }
//+ }
//+ void	displayVariables(Domain *domain)
//+ {	int
//+ 		ivar=domain->getDisplayVariable(),
//+ 		icomp=domain->getDisplayVariableComp();
//+ 	class Variable	*var;
//+ 	if (ivar>=0)
//+ 	{	var=domain->variable+ivar;
//+ 		if (var->type==pores)
//+ 			glPointSize(P.disp[size]);
//+ 		else
//+ 			glPointSize(scene.node[size]);
//+ 		glDisable(GL_LIGHTING);
//+ 		if(icomp<0) showVectors(ivar,domain);
//+ 		else
//+ 		{	double	vmn,vmx;
//+ 			if (colorscale.relative)
//+ 			{	var->computeLimits();
//+ 				var->getLimits(icomp,vmn,vmx);
//+ 				if(var->type==pores)
//+ 				{	P.disp[vmin]=vmn;
//+ 					P.disp[vmax]=vmx;
//+ 				}
//+ 				else
//+ 				{	scene.node[vmin]=vmn;
//+ 					scene.node[vmax]=vmx;
//+ 				}
//+ 			}
//+ 			else
//+ 			{	if(var->type==pores)
//+ 				{	vmn=P.disp[vmin];
//+ 					vmx=P.disp[vmax];
//+ 				}
//+ 				else
//+ 				{	vmn=scene.node[vmin];
//+ 					vmx=scene.node[vmax];
//+ 				}
//+ 			}
//+ 			if(vmx>vmn)
//+ 				showScalar
//+ 				(	ivar,icomp,
//+ 					vmn,vmx,
//+ 					domain
//+ 				);
//+ 			else
//+ 			{	if (var->dimension>1)
//+ 				fprintf
//+ 				(	stderr,
//+ 					"Component %d of ",icomp
//+ 				);
//+ 				if (vmx==vmn)
//+ 				fprintf
//+ 				(	stderr,
//+ 					"%s has a constant value of %g in the whole domain %s\n",
//+ 					var->name,vmn,domain->name
//+ 				);
//+ 				else
//+ 				{	fprintf
//+ 					(	stderr,
//+ 						"%s limits are in a wrong range: min=%g, max=%g\n",
//+ 						var->name,scene.node[vmin],scene.node[vmax]
//+ 					);BUG("Internal inconsistency");
//+ 				}
//+ 			}
//+ 		}
//+ 		glEnable(GL_LIGHTING);
//+ 	}
//+ }
void	displayAxes() {
	GLfloat mat[4] = { 0.6, 0.8, 0.3, 1.0 },
		shine=0.6;
	mat[0]=.8;mat[1]=0.7;mat[2]=0.5;
	glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
	glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
	glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
	if (!showpar[showAxes]) return;
	glEnable(GL_LIGHTING);
	glLineWidth (axes[axeswidth]);
	glBegin (GL_LINES);
	glVertex3f(xo[0],xo[1],xo[2]);
	glVertex3f(xo[0]+axes[xaxislength]*lx,xo[1],xo[2]);
	glVertex3f(xo[0],xo[1],xo[2]);
	glVertex3f(xo[0],xo[1]+axes[yaxislength]*ly,xo[2]);
	glVertex3f(xo[0],xo[1],xo[2]);
	glVertex3f(xo[0],xo[1],xo[2]+axes[zaxislength]*lz);
	glEnd ();
	
	//Added
	glPushMatrix();
	
	glTranslatef (xo[0], xo[1], xo[2]+axes[zaxislength]*lz);
	glutSolidCone (axes[arrowidth]*lz, axes[arrowheight]*lz, 4, 4);
	//glTranslatef (-xo[0], -xo[1], -(xo[2]+axes[zaxislength]*lz));
	
	glPopMatrix();
}
/* ARGSUSED1 */
void	getScaling
(
	double	&l0,
	double	&l1
)
{
	l0=lmin,l1=lmax;
}
void	getXLimits
(	double	*x0,
	double	*x1
)
{
	for (int i=0; i<3; i++)
	{	x0[i]=xmin[i];
		x1[i]=xmax[i];
	}
}
void	Exit()
{
	char	s[MAXLINLEN];
	if (animation) {
		animation=0;
	}
	printf("Exit (y/n)? ");fflush(stdout);
	fgets(s,MAXLINLEN,stdin);
	if (*s!='y') return;
//-	cleanup();
	exit(0);
}
void	Quit()
{
	char	s[MAXLINLEN];
	if (animation) {
		animation=0;
	}
	printf("Quit (y/n)? ");fflush(stdout);
	fgets(s,MAXLINLEN,stdin);
	if (*s!='y') return;
//-	cleanup();
	exit(0);
}
void	menu(int value) {
	printf("Menu %d: ",value);
	switch (value) {
		case 0:
			showpar[showBoundaryVertexes] = showpar[showBoundaryVertexes]==0?1:0;
			printf("Show Boundary Vertexes = %d\n",showpar[showBoundaryVertexes]);
			break;
		case 1:
			showpar[showBoundaryGrid] = showpar[showBoundaryGrid]==0?1:0;
			printf("Show Boundary Grid = %d\n",showpar[showBoundaryGrid]);
			break;
		case 2:
			showpar[showBoundaryFaceCenters] = showpar[showBoundaryFaceCenters]==0?1:0;
			printf("Show Boundary Faces = %d\n",showpar[showBoundaryFaceCenters]);
			break;
		case 3:
			showpar[showNodes] = showpar[showNodes]==0?1:0;
			printf("Show Grid Nodes = %d\n",showpar[showNodes]);
			break;
		case 4:
			showpar[showFaceCenters] = showpar[showFaceCenters]==0?1:0;
			printf("Show Face Centers = %d\n",showpar[showFaceCenters]);
			break;
		case 5:
			showpar[showCellCenters] = showpar[showCellCenters]==0?1:0;
			printf("Show Cell Centers = %d\n",showpar[showCellCenters]);
			break;
		case 6:
			showpar[showGrid] = showpar[showGrid]==0?1:0;
			printf("Show Grid = %d\n",showpar[showGrid]);
			break;
		case 7:
			printf("Option 7 Not implemented \n");
			break;
		case 8:
			showpar[showBonds] = showpar[showBonds]==0?1:0;
			printf("Show Bonds = %d\n",showpar[showBonds]);
			break;
		case 9:
			showpar[showVariables] = showpar[showVariables]==0?1:0;
			printf("Show Variables = %d\n",showpar[showVariables]);
			break;
		case 10:
			showpar[showAxes] = showpar[showAxes]==0?1:0;
			printf("Show Axes = %d\n",showpar[showAxes]);
			break;
		case 11:
			showpar[showBoundaryVectors] = showpar[showBoundaryVectors]==0?1:0;
			printf("Show Boundary Vectors = %d\n",showpar[showBoundaryVectors]);
			break;
		case 12:
//-			selectVariable();
			break;
		case 13:
//-			init();
			break;
		case 14:
			if (finished) break;
			animation=animation==0?1:0;
			if (animation) 
				glutIdleFunc(animate);
			else
			{	glutIdleFunc(NULL);
				printf("\n");fflush(stdout);
			}
		break;
		case 15:
			helpDisplay();
		break;
		case 16:
		printf(ABOUT);
		break;
		case 17:
		printf("\nSend your comments to andrei@smirnov.mae.wvu.edu\n");
//-		cleanup();
		exit(0);
		break;
	}
	fflush(stdout);
	//glutPostRedisplay();
}
void displayMenu()
{
	printf("----------------- MENU -----------------\n");
	printf(" ACTION                   (key-shortcut)\n");
	printf(" 0. Quit menu             (<Enter>)\n");
	printf("Toggle the display of        \n");
	printf(" 1. Boundary vertexes        \n");
	printf(" 2. Boundary grid         (g)\n");
	printf(" 3. Boundary faces           \n");
	printf(" 4. Nodes                 (n)\n");
	printf(" 5. Frame                 (f)\n");
//	printf(" 6. Cells                 (c)\n");
//	printf(" 7. Edges                 (e)\n");
//	printf(" 8. Bonds                 (b)\n");
//	printf(" 9. Variables             (V)\n");
	printf("10. Axes                  (a)\n");
	printf("11. Boundary Vectors      (v)\n");
	printf("Actions:\n");
//	printf("12.                          \n");
//	printf("13.                          \n");
	printf("14. Initialize            (r)\n");
	printf("15. Run/Stop              (r)\n");
	printf("16. Help                  (?)\n");
	printf("17. About                    \n");
	printf("18. Quit program          (q)\n");
}
void consoleMenu()
{	int	selection=-1;
	do
	{	char	buf[MAXLINLEN];
		displayMenu();
		printf("Select: ");
		selection=0;
		fgets(buf,MAXLINLEN,stdin);
		sscanf(buf,"%d",&selection);
		if (selection==0) return; 
		if (selection>0&&selection<=15)
			menu(selection-1);
		else printf("Wrong selection: %d\n",selection);
		//else	break;
	}	while (selection!=0);
}
void	setBackgroundRun()
{
	showpar[showRun]=0;
}
void	setForegroundRun()
{
	showpar[showRun]=1;
}
void	toggleSpheres() {
	showpar[showSpheres]=showpar[showSpheres]?0:1;
	printf("Show Spheres = %d\n",showpar[showSpheres]);
}
void	switchColorScheme() {
	scene.color.scheme=(ColorSchemes)(((int)scene.color.scheme+1)%(int)maxColorSchemes);
	printf("Color Scheme = %d\n",(int)scene.color.scheme);
}
void	runmany()
{	int n=0;
	printf("Enter number or iterations: ");fflush(stdout);
	scanf("%d",&n);
	printf("Executing %d iterations ...",n);fflush(stdout);
	domain->run(n);
	printf(" done\n");
}
//-void	setWriteFormat()
//-{	using namespace Output;
//-	int	selection=-1;
//-	do
//-	{	char	buf[MAXLINLEN];
//-		printf("Output formats:\n");
//-		for(int i=0;i<maxoutypes;i++)
//-			printf("%d - %s\n",i,Output::outypename[i]);
//-		printf("Select: ");
//-		selection=0;
//-		fgets(buf,MAXLINLEN,stdin);
//-		sscanf(buf,"%d",&selection);
//-		if (selection>=0&&selection<maxoutypes)break;
//-		printf("Wrong selection: %d\n",selection);
//-	}	while (selection!=0);
//-	outype=(OutputTypes)selection;
//-	printf("Selected output type: %s\n",Output::outypename[outype]);fflush(stdout);
//-}
//-void	writeData() {
//-	char	*s,buf[MAXLINLEN+1],
//-		outfilename[MAXLINLEN];
//-	using namespace	Output;
//-	sprintf(outfilename,"%s",domain->name);
//-	printf("Write:\n\t grid (g)\n\t boundary (b)\n\t current-variable (v)\n\t field-data (f)\n\t pore-data (p)?\nEnter your choice: ");fflush(stdout);
//-	fgets(buf,MAXLINLEN,stdin);
//-	for (s=buf; isspace(*s); s++);
//-	switch (*s)
//-	{
//-		case 'g':
//-		Domain::Save(outfilename);
//-		break;
//-		case 'b':
//-		printf("Case b: outputBoundary: NOT IMPLEMENTED\n");
//-		break;
//-		case 'v':
//-		printf("Case v: NOT IMPLEMENTED\n");
//-		break;
//-		case 'f':
//-		printf("NOT IMPLEMENTED\n");
//-		break;
//-		case 'p':
//-		printf("Ouput pores not implemented yet\n");
//-		break;
//-		default:
//-		printf("Operation skipped");
//-	}
//-	printf("\n");
//-}
void	toggleWindowDump() {
	printf("Enter window dump time interval: ");fflush(stdout);
	scanf("%g",&wdtime);
	if(wdtime>0.0)
	{	showpar[dumpWindow]=1;
		iwdump=0;nwdump=(int)(wdtime/Run::time.step+0.5);
		printf("Window dump = %d, dump time interval = %g\n",showpar[dumpWindow],wdtime);
	}
	else
	{	showpar[dumpWindow]=0;
		printf("Window dump disabled\n",showpar[dumpWindow],wdtime);
	}
}
void	dumpwindow()
{	using namespace IO;
	xwd(windowname);
}
void	commandMode()
{
	static struct Command
	{
		char	*name;
		void	(*f)();
	}	command[]	= 
	{
//-		{"v",selectVariable},
		{(char*)"m",consoleMenu},
//-		{"w",writeData},
//-		{"wf",setWriteFormat},
//-		{"lim",printGridXLimits},
//-		{"gvec",printGridVecLimits},
		{(char*)"bg",setBackgroundRun},
		{(char*)"cs",switchColorScheme},
		{(char*)"fg",setForegroundRun},
		{(char*)"?",helpCommand},
		{(char*)"h",helpCommand},
		{(char*)"help",helpCommand},
		{(char*)"r",runmany},
		{(char*)"sp",toggleSpheres},
		{(char*)"wd",toggleWindowDump},
		{(char*)"dw",dumpwindow},
		{(char*)"e",Exit},
		{(char*)"q",Quit},
		{(char*)".",NULL},
		{(char*)"\0",NULL}
	},	*cmd;
	char	buf[MAXLINLEN+1],*p,*s;
	printf("Command mode\n");
	while (1) {
		printf(":");fflush(stdout);
#ifdef WITH_DOVE
	if(iproc==0)
#endif
	{ fgets(buf,MAXLINLEN,stdin);
	}
/*
#ifdef WITH_MPI
		COMM_WORLD.Bcast(buf, MAXLINLEN, MPI_CHAR, 0);
#endif
*/
		for (s=buf; isspace(*s); s++);
		if ((p=strchr(s,'\n'))!=NULL) *p='\0';
		for (cmd=command; *cmd->name!='\0'; cmd++)
			if (strcmp(s,cmd->name)==0) break;
		switch (*cmd->name)
		{
			case '\0':
			if (*s=='\0')
			{	printf("Display-mode\n");fflush(stdout);
				return;
			}
			else
				printf("Unknown command\n");
			break;
			case '.':
				printf("\nDisplay-mode\n");fflush(stdout);
				return;
			default:
			(*cmd->f)();
		}
	}
}
void	display() {
	//INITIALIZE GL:
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	// added
	// Apply the camera translations
	glTranslatef(0.0, 0.0, -zoom);
	glTranslatef(tx, ty, 0.0);
	glRotatef(rotx, 1.0, 0.0, 0.0);
	glRotatef(roty, 0.0, 1.0, 0.0);
	
	displayAxes();
	//	DISPLAY DOMAIN:
	refresh();
	if (showpar[showGrid]) {
		cout << "showGrid disabled\n";
//-		REAL	color[3];
//-		//int showboundarygrid=showpar[showBoundaryGrid];
//-		glDisable(GL_LIGHTING);
//-		glLineWidth (scene.line[thickness]);
//-		getElementColor(edges,color);
//-		glColor3f(color[red], color[green], color[blue]);
//-		List<Cell> *cells=domain->Cells();
//-		if (cells->number()>0) {
//-			const int dim1=DIM+1;
//-			if(Run::option.debug)cout<<"Displaying "<<cells->number()<<" cells\n";
//-			int icell=0;
//-			cells->goFirst();
//-			do {
//-				if(Run::option.debug){cout<<icell++<<endl;cout.flush();}
//-				Cell *cell=cells->Current();
//-				Node **verts=cell->Vertexes();
//-				glBegin(GL_LINES);
//-				for(int step=1;step<3;step++)//! for DIM=4 only
//-				for(int iv=0;iv<dim1;iv++) {
//-					Node 
//-						*v0=verts[iv],
//-						*v1=verts[(iv+step)%dim1]; 
//-					REAL 
//-						*x=v0->Coordinates(),
//-						*y=v1->Coordinates();
//-					if(Run::option.debug)
//-						cout<<iv<<" -> "<<(iv+step)%dim1
//-						<<":  "<<x[0]<<' '<<x[1]<<' '<<x[2]
//-						<<"-> "<<y[0]<<' '<<y[1]<<' '<<y[2]<<endl;
//-					glVertex3f(x[0],x[1],x[2]);
//-					glVertex3f(y[0],y[1],y[2]);
//-				}
//-				glEnd();
//-				cells->goNext();
//-			}	while(!cells->isFirst());
//-		}
//-		glEnable(GL_LIGHTING);
	}
	if (showpar[showBoundaryGrid]) {
		cout << "showBoundaryGrid disabled\n";
//+		REAL	color[3];
//+		glDisable(GL_LIGHTING);
//+		glLineWidth (scene.mesh.line[thickness]);
//+		getElementColor(boundary_edges,color);
//+		glColor3f(color[red], color[green], color[blue]);
//+		REAL xm[2][DIM], edges[12][2][3]; 
//+		for (int i=0;i<DIM;i++) {
//+			xm[0][i]=domain->minBound(i);
//+			xm[1][i]=domain->maxBound(i);
//+		}
//+		//0: 0 0 0 - 0 0 1
//+		edges[0][0][0]=xm[0][0];
//+		edges[0][0][1]=xm[0][1];
//+		edges[0][0][2]=xm[0][2];
//+		edges[0][1][0]=xm[0][0];
//+		edges[0][1][1]=xm[0][1];
//+		edges[0][1][2]=xm[1][2];
//+		//1: 0 0 0 - 0 1 0
//+		edges[1][0][0]=xm[0][0];
//+		edges[1][0][1]=xm[0][1];
//+		edges[1][0][2]=xm[0][2];
//+		edges[1][1][0]=xm[0][0];
//+		edges[1][1][1]=xm[1][1];
//+		edges[1][1][2]=xm[0][2];
//+		//2: 0 0 0 - 1 0 0
//+		edges[2][0][0]=xm[1][0];
//+		edges[2][0][1]=xm[0][1];
//+		edges[2][0][2]=xm[0][2];
//+		edges[2][1][0]=xm[0][0];
//+		edges[2][1][1]=xm[0][1];
//+		edges[2][1][2]=xm[0][2];
//+		// 
//+		//3: 0 0 1 - 0 1 1
//+		edges[3][0][0]=xm[0][0];
//+		edges[3][0][1]=xm[0][1];
//+		edges[3][0][2]=xm[1][2];
//+		edges[3][1][0]=xm[0][0];
//+		edges[3][1][1]=xm[1][1];
//+		edges[3][1][2]=xm[1][2];
//+		//4: 0 0 1 - 1 0 1
//+		edges[4][0][0]=xm[0][0];
//+		edges[4][0][1]=xm[0][1];
//+		edges[4][0][2]=xm[1][2];
//+		edges[4][1][0]=xm[1][0];
//+		edges[4][1][1]=xm[0][1];
//+		edges[4][1][2]=xm[1][2];
//+		// 
//+		//5: 0 1 0 - 1 1 0
//+		edges[5][0][0]=xm[0][0];
//+		edges[5][0][1]=xm[1][1];
//+		edges[5][0][2]=xm[0][2];
//+		edges[5][1][0]=xm[1][0];
//+		edges[5][1][1]=xm[1][1];
//+		edges[5][1][2]=xm[0][2];
//+		//6: 0 1 0 - 0 1 1
//+		edges[6][0][0]=xm[0][0];
//+		edges[6][0][1]=xm[1][1];
//+		edges[6][0][2]=xm[0][2];
//+		edges[6][1][0]=xm[0][0];
//+		edges[6][1][1]=xm[1][1];
//+		edges[6][1][2]=xm[1][2];
//+		// 
//+		//7: 1 0 0 - 1 0 1
//+		edges[7][0][0]=xm[1][0];
//+		edges[7][0][1]=xm[0][1];
//+		edges[7][0][2]=xm[0][2];
//+		edges[7][1][0]=xm[1][0];
//+		edges[7][1][1]=xm[0][1];
//+		edges[7][1][2]=xm[1][2];
//+		//8: 1 0 0 - 1 1 0
//+		edges[8][0][0]=xm[1][0];
//+		edges[8][0][1]=xm[0][1];
//+		edges[8][0][2]=xm[0][2];
//+		edges[8][1][0]=xm[1][0];
//+		edges[8][1][1]=xm[1][1];
//+		edges[8][1][2]=xm[0][2];
//+		// 
//+		//9: 0 1 1 - 1 1 1
//+		edges[9][0][0]=xm[0][0];
//+		edges[9][0][1]=xm[1][1];
//+		edges[9][0][2]=xm[1][2];
//+		edges[9][1][0]=xm[1][0];
//+		edges[9][1][1]=xm[1][1];
//+		edges[9][1][2]=xm[1][2];
//+		//10: 1 1 0 - 1 1 1
//+		edges[10][0][0]=xm[1][0];
//+		edges[10][0][1]=xm[1][1];
//+		edges[10][0][2]=xm[0][2];
//+		edges[10][1][0]=xm[1][0];
//+		edges[10][1][1]=xm[1][1];
//+		edges[10][1][2]=xm[1][2];
//+		// 
//+		//11: 1 1 1 - 1 0 1
//+		edges[11][0][0]=xm[1][0];
//+		edges[11][0][1]=xm[1][1];
//+		edges[11][0][2]=xm[1][2];
//+		edges[11][1][0]=xm[1][0];
//+		edges[11][1][1]=xm[0][1];
//+		edges[11][1][2]=xm[1][2];
//+
//+		glBegin(GL_LINES);
//+		for(int ie=0;ie<12;ie++) {
//+			for(int ivert=0;ivert<2;ivert++) {
//+				REAL *x=edges[ie][ivert];
//+				glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//+			}//end for(ivert)
//+		}//end for(ie)
//+		glEnd();
	}
	if(showpar[showFrame]) {
		REAL	edgescolor[3];
		getElementColor(frame,edgescolor);
		glDisable(GL_LIGHTING);
		glLineWidth (scene.frame.line[thickness]);
		glColor3f(edgescolor[red], edgescolor[green], edgescolor[blue]);
		REAL v[2][2][2][3];
		//8 vertexes of a box:
		v[0][0][0][0]=xmin[0];
		v[0][0][0][1]=xmin[1];
		v[0][0][0][2]=xmin[2];

		v[0][0][1][0]=xmin[0];
		v[0][0][1][1]=xmin[1];
		v[0][0][1][2]=xmax[2];

		v[0][1][0][0]=xmin[0];
		v[0][1][0][1]=xmax[1];
		v[0][1][0][2]=xmin[2];

		v[0][1][1][0]=xmin[0];
		v[0][1][1][1]=xmax[1];
		v[0][1][1][2]=xmax[2];

		v[1][0][0][0]=xmax[0];
		v[1][0][0][1]=xmin[1];
		v[1][0][0][2]=xmin[2];

		v[1][0][1][0]=xmax[0];
		v[1][0][1][1]=xmin[1];
		v[1][0][1][2]=xmax[2];

		v[1][1][0][0]=xmax[0];
		v[1][1][0][1]=xmax[1];
		v[1][1][0][2]=xmin[2];

		v[1][1][1][0]=xmax[0];
		v[1][1][1][1]=xmax[1];
		v[1][1][1][2]=xmax[2];

		glBegin(GL_LINES);
		for(int i=0;i<2;i++)
		for(int j=0;j<2;j++) {
			glVertex3f((GLfloat)v[i][j][0][0],(GLfloat)v[i][j][0][1],(GLfloat)v[i][j][0][2]);
			glVertex3f((GLfloat)v[i][j][0][0],(GLfloat)v[i][j][0][1],(GLfloat)v[i][j][1][2]);
			glVertex3f((GLfloat)v[i][0][j][0],(GLfloat)v[i][0][j][1],(GLfloat)v[i][0][j][2]);
			glVertex3f((GLfloat)v[i][0][j][0],(GLfloat)v[i][1][j][1],(GLfloat)v[i][0][j][2]);
			glVertex3f((GLfloat)v[0][i][j][0],(GLfloat)v[0][i][j][1],(GLfloat)v[0][i][j][2]);
			glVertex3f((GLfloat)v[1][i][j][0],(GLfloat)v[0][i][j][1],(GLfloat)v[0][i][j][2]);
		}
		glEnd();
	}
	if(showpar[showNodes]) {
		REAL	ncolor[3],bcolor[3];
		int showboundarynodes=showpar[showBoundaryVertexes];
		getElementColor(faces,bcolor);
		if(scene.color.scheme==colorByTime)Palette::init(Run::time.start,Run::time.end);
		if(scene.color.scheme==colorByType)Palette::init(scene.color.minvalue,scene.color.maxvalue);
		if(scene.color.scheme==colorByMass)Palette::init(scene.color.minvalue,scene.color.maxvalue);
		getElementColor(nodes,ncolor);
		if(showpar[showSpheres]) {// Show as spheres
			if(Run::option.debug)printf("Show Nodes as Spheres\n");
			GLfloat	radius=scene.mesh.node[size];
			if(radius<=0.0)radius=0.5*Potential::lengthscale();
			GLfloat mat[4],//- = { 0.6, 0.8, 0.3, 1.0 },
				shine=0.6;
			mat[3]=1.0;
			glEnable(GL_LIGHTING);
			glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
			if (domain->Molecules()->number()>0) {
//?				List<Molecule> *nodes=domain->Molecules();
				Collection<Molecule> *nodes=domain->Molecules();
				nodes->goFirst();
				do {
					Molecule *node=nodes->Current();
					REAL
						size=species[node->Type()].Size(),
						*x=node->Coordinates();
					if(size>0.0) radius=0.5*size;
					glTranslatef ((GLfloat)x[0], (GLfloat)x[1], (GLfloat)x[2]); 
					switch(scene.color.scheme) {
						case colorByBoundary:
//-						if (node->atBoundary()==TRUE) {
//-							mat[0]=bcolor[red];mat[1]=bcolor[green];mat[2]=bcolor[blue];
//-						} else {
							mat[0]=ncolor[red];mat[1]=ncolor[green];mat[2]=ncolor[blue];
//-						}
						break;
						case colorByType:
							Palette::pickcolor((REAL)node->Type(),ncolor);
							mat[0]=ncolor[red];mat[1]=ncolor[green];mat[2]=ncolor[blue];
						break;
						case colorByTime:
							Palette::pickcolor(Run::time.current,ncolor);
							mat[0]=ncolor[red];mat[1]=ncolor[green];mat[2]=ncolor[blue];
						break;
						case colorByMass:
							Palette::pickcolor(species[node->Type()].Mass(),ncolor);
							mat[0]=ncolor[red];mat[1]=ncolor[green];mat[2]=ncolor[blue];
						break;
						default:
							mat[0]=ncolor[red];mat[1]=ncolor[green];mat[2]=ncolor[blue];
						break;
					}
					mat[3]=(x[3]-xmin[3])/(xmax[3]-xmin[3]);//set by time
					glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
					glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
					glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
					glutSolidSphere (radius, 12, 12);
					glTranslatef (-(GLfloat)x[0],-(GLfloat)x[1],-(GLfloat)x[2]); 
					nodes->goNext();
				}	while(!nodes->isFirst());
			}
		} else {// Show nodes as points:
			if(Run::option.debug)printf("Show Nodes as Points\n");
			glPointSize(scene.mesh.node[size]);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			if (domain->Molecules()->number()>0) {
//?				List<Molecule> *nodes=domain->Molecules();
				Collection<Molecule> *nodes=domain->Molecules();
				nodes->goFirst();
				do {
					Molecule *node=nodes->Current();
					REAL *x=node->Coordinates();
					switch(scene.color.scheme) {
						case colorByBoundary:
//-						if (node->atBoundary()==TRUE) {
//-							glColor3f(bcolor[red], bcolor[green], bcolor[blue]);
//-						} else {
//-						}
						break;
						case colorByTime:
							Palette::pickcolor(Run::time.current,ncolor);
						break;
						case colorByType:
							Palette::pickcolor((REAL)node->Type(),ncolor);
						break;
						case colorByMass:
							Palette::pickcolor(species[node->Type()].Mass(),ncolor);
						break;
					}
					glColor3f(ncolor[red], ncolor[green], ncolor[blue]);
					glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
					nodes->goNext();
				}	while(!nodes->isFirst());
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
	glPopMatrix();
//	glFlush();
//glutSwapBuffers();
	glXSwapBuffers(dpy, win);
}
void	animate() {
	if (domain->run(1)==0) {
//- finish();
		finished=1;
		animation=0;
///DDD		glutIdleFunc(NULL);//Causes segmentation fault (!?)
		display();
		return;
	}
	if(!Run::option.xterm)printf("\rTime = %g%c[K",Run::time.current,ESC);fflush(stdout);
	if (showpar[showRun])
	{///	refresh();
		display();
	}
	if (showpar[dumpWindow]) 
	{
		if (iwdump++%nwdump==0)
			dumpwindow();
	}
}
void	helpCommand() {
	printf("e     - save and exit\n");
	printf("h     - show help\n");
	printf("bg    - set background run (no automatic redisplay)\n");
	printf("dw    - dump window\n");
	printf("fg    - set foreground run (automatic redisplay)\n");
	printf("lim   - display domain limits\n");
	printf("m     - select from a menu\n");
	printf("gvec  - display grid vector limits\n");
	printf("gv    - display grid variables limits\n");
	printf("pv    - display pore variables limits\n");
	printf("q     - quit without saving\n");
	printf("r     - run\n");
	printf("u     - uniform\n");
	printf("var   - select a variable\n");
	printf("w     - write data at this time step into %s-<number>.dat.gz\n",Run::outputname);
	printf("wd    - toggle window dump\n");
	printf(".     - quit command mode\n");
}
void	reshape(int w, int h)
{
	double
		lmin,lmax;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	getScaling(lmin,lmax);
	gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 0.1*lmin, 10*lmax);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	/* REAL	
		r0=5.,
		r1=10.0;
	double
		lmin,lmax,
		xmin[DIM],xmax[DIM];
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w <= h) 
	getXLimits(xmin,xmax);
	glOrtho 
	(	xmin[0], 
		xmax[0], 
		xmin[1], 
		xmax[1], 
		xmin[2], 
		xmax[2]
	);
	glLoadIdentity ();
	getScaling(lmin,lmax);
//	gluPerspective (45., (GLdouble)w/(GLdouble)h, .1*lmin, 10*lmax);
  	gluPerspective (50., (GLdouble)w/(GLdouble)h, 0.1*lmin, 10*lmax);
	glMatrixMode(GL_MODELVIEW); */
}
void	mouse(int button, int state, int x, int y)
{
	lastx = x;
	lasty = y;
	
	mouseButtons[0] = 0;
	mouseButtons[1] = 0;
	mouseButtons[2] = 0;
	
	switch(button)
	{
		case LEFT_BUTTON:
		{
			mouseButtons[0] = 1;
		}break;
		
		case MIDDLE_BUTTON:
		{
			mouseButtons[1] = 1;
		}break;
		
		case RIGHT_BUTTON:
		{
			mouseButtons[2] = 1;
		}break;
	}
}
void	motion(int x, int y) {
	
	REAL diffx = x - lastx;
	REAL diffy = y - lasty;
	
	lastx = x;
	lasty = y;
	
	if(mouseButtons[0])
	{
		rotx += 0.09f * diffy;
		roty += 0.09f * diffx;
	}
	else if(mouseButtons[1])
	{
		zoom += 0.05f * diffy;
	}
	else if(mouseButtons[2])
	{
		tx += 0.05f * diffx;
		ty -= 0.05f * diffy;
	}
}
/* ARGSUSED3 */
void	keyboard(unsigned int key) {
	using namespace IO;
	if (Run::option.debug) printf("key=%c (%d)\n",(char)key,key);
	if (key<32||key>127)return;
	switch (key) 
	{
		case '+':
			if (finished) break;
			if (domain->run(1)==0)	finished=1;
			printf("Time: %-10.3g\n",Run::time.current);fflush(stdout);
			break;
		case '-':
			if (finished) break;
			if (domain->run(-1)==0)	finished=1;
			printf("Time: %-10.3g\n",Run::time.current);fflush(stdout);
			break;
		case 'A':
		case 'a':
			showpar[showAxes] = showpar[showAxes]==0?1:0;
			printf("Show Axes = %d\n",showpar[showAxes]);
			break;
		case 'B':
			showpar[showBoundaryFaces] = showpar[showBoundaryFaces]==0?1:0;
			printf("Show Boundary Faces = %d\n",showpar[showBoundaryFaces]);
			break;
		case 'c':
		case 'C':
			switchColorScheme();
			break;
//-			domain->domain2D();
//-			printf("\r2D domaining = %d nodes\n",domain->Nodes()->number());fflush(stdout);
//-			printf("\rEnergy = %g      ",domain->Energy());fflush(stdout);
//-		break;
//-		case 'C':
//-			domain->domain3D();
//-			printf("\r2D domaining = %d nodes\n",domain->Nodes()->number());fflush(stdout);
//-			printf("\rEnergy = %g      ",domain->Energy());fflush(stdout);
//-		break;
		case 'G':
			showpar[showGrid] = showpar[showGrid]==0?1:0;
			printf("Show Grid = %d\n",showpar[showGrid]);
			break;
		case 'g':
			showpar[showBoundaryGrid] = showpar[showBoundaryGrid]==0?1:0;
			printf("Show Boundary Grid = %d\n",showpar[showBoundaryGrid]);
			break;
		case 'f':
			showpar[showFrame] = showpar[showFrame]==0?1:0;
			printf("Show Frame = %d\n",showpar[showFrame]);
			break;
		case 'I':
		case 'i':
//-			init();
			break;
		case 'M':
		case 'm':
			consoleMenu();
			break;
		case 'N':
			showpar[showNodes] = showpar[showNodes]==0?1:0;
			printf("Show Grid Nodes = %d\n",showpar[showNodes]);
			//glutPostRedisplay();
			break;
		case 'n':
			showpar[showBoundaryVertexes] = showpar[showBoundaryVertexes]==0?1:0;
			printf("Show Boundary Vertexes = %d\n",showpar[showBoundaryVertexes]);
			break;
		case 'l':
		case 'L':
			readconf();
			break;
		case 'v':
			showpar[showBoundaryVectors] = showpar[showBoundaryVectors]==0?1:0;
			printf("Show Boundary Vectors = %d\n",showpar[showBoundaryVectors]);
			break;
		case 'R':
		case 'r':
			if (finished) break;
			animation=animation==0?1:0;
			if (Run::option.verbose|Run::option.debug)
				printf("animation=%d\n",animation);
			if(Run::option.xterm)printf("%c[2J",ESC);//erase screen
			break;
		case 'S':
			finished=0;
		case 's':
			if (finished) break;
			if(Run::option.xterm)printf("%c[2J",ESC);//erase screen
			if (domain->run(1)==0) {
//-	finish();
				finished=1;
			}
			///refresh();
			printf("Time: %-9.3g\n",Run::time.current);fflush(stdout);
			break;
//-		case 't':
//-			domain->tagBoundaryNodes();
//-			printf("Time: %g        \n",Run::time.current);fflush(stdout);
//-		break;
//-		case 'u':
//-			domain->relax();
//-			printf("Time: %g        \n",Run::time.current);fflush(stdout);
//-		break;
		case 'w':
		//-	dumpwindow();
		{	static int idump=0;
			char filename[MAXLINLEN];
			sprintf(filename,"%s_dump%d",Run::outputname,idump++);
			domain->save(filename);
		}
		break;
		case 'z':
			zoom -= 1.5;
			break;
		case 'Z':
			zoom += 1.5;
			break;
		case ':':
		case '.':
			commandMode();
			break;
		case 27:
			Exit();
			break;
		case '?':
			helpDisplay();
			break;
	}
}
static int
events_loop(Display *dpy, Window win) 
{	XEvent event;
	do 
	{	char buf[31];
		KeySym keysym;
		XNextEvent(dpy, &event);
		switch(event.type) 
		{	case Expose:
				break;
			case ConfigureNotify: 
				reshape(event.xconfigure.width, event.xconfigure.height);
			//	{	/* this approach preserves a 1:1 viewport aspect ratio */
			//		int vx, vy, vw, vh;
			//		int ew = event.xconfigure.width, eh = event.xconfigure.height;
			//		if (ew >= eh) 
			//		{
			//			vx = 0;
			//			vy = (eh - ew) >> 1;
			//			vw = vh = ew;
			//		} 
			//		else 
			//		{
			//			vx = (ew - eh) >> 1;
			//			vy = 0;
			//			vw = vh = eh;
			//		}
			//		glviewport(vx, vy, vw, vh);
			//	}
				break;
			case KeyPress:
				(void) XLookupString(&event.xkey, buf, sizeof(buf), &keysym, NULL);
				keyboard((unsigned char)keysym);
				break;
			case ButtonPress:
				mouse(event.xbutton.button,GLUT_DOWN,event.xbutton.x,event.xbutton.y);
				break;
			case MotionNotify:
				motion(event.xbutton.x, event.xbutton.y);
				break;
			default:
				break;
		}
		if (animation)
		do
		{	animate();
			//display();
		}	while(!XPending(dpy));
	} while (XPending(dpy));
	display();
	return 0;
}
void	run() {
	while (events_loop(dpy, win)==0);
}
void displaymessage(char *msg)
{
	puts(msg);
}
void	setElementColor(int element_type, REAL rgbcolor[])
{
	if (element_type<maxelements)
	for (int irgb=0; irgb<3; irgb++)
		disp[element_type].rgbcolor[irgb]=rgbcolor[irgb];
}
void	getElementColor(int element_type, REAL rgbcolor[])
{
	if (element_type<maxelements)
		for (int irgb=0; irgb<3; irgb++)
			rgbcolor[irgb]=disp[element_type].rgbcolor[irgb];
	else {
		fprintf(stderr,"getElementColor: no color for element type %d\n",element_type);
		fflush(stderr);
	}
}
} //END NAMESPACE Gui
namespace Palette {
	REAL
		vmn,
		vmx,
		vav,
		dvi;
	void init(REAL vmin, REAL vmax) {
		vmn=vmin;
		vmx=vmax;
		vav=.5*(vmn+vmx);
		dvi=2.0/(vmx-vmn);
	//-				alow =PI/(2.*(vav-vmn));
	//-				ahigh=PI/(2.*(vmx-vav));
	}
	void pickcolor(REAL var, REAL color[]) {
		double	r,g,b;
	//-			double	angle;
		if (var<=vmn)
		{	r=g=0.0; b=1.0; }
		else
		if (var>=vmx)
		{	r=1.0; g=b=0.0; }
		if (var<vav) {	
	//-			angle=(double)(var-vmn)*alow;
	//-			b=cos(angle); g=sin(angle);
			g=(var-vmn)*dvi;
			b=1.0-g;
			r=0.0;
		} else {	
	//-			angle=(double)(var-vav)*ahigh;
	//-			g=cos(angle); r=sin(angle);
			r=(var-vav)*dvi;
			g=1.0-r;
			b=0.0;
		}
		color[0]=(REAL)r;
		color[1]=(REAL)g;
		color[2]=(REAL)b;
	}
}// END NAMESPACE Palette

