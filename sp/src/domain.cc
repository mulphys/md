/*! \file domain.cc
 * \brief Implementation of Domain and Boundary classes
 */
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <limits>
using namespace std;
#include "def.h"
#include "io.h"
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


void Boundary::init(int i, int j, REAL ymin[], REAL ymax[], int nspecies) {
	idir=i; iside=j; xmin=ymin; xmax=ymax;
	reactions=new Reaction[nspecies];
}
//! Injecting a specie from accross the boundary with a specified velocity
void Boundary::Inject(
	int ispecie, 
	REAL vavx,//average velocity in on direction
	REAL velx,//actual velocity normal to the boundary (randomized)
	REAL temperature,
	Collection<Molecule> *molecules
) {
	Molecule *molecule=molecules->append();
	molecule->Type(ispecie);
	//Select position: 
	if(iside==0) {
		molecule->Coordinate(idir,xmin[idir]+numeric_limits<REAL>::epsilon());
//-		molecule->Velocity(idir,RND*velx);
		molecule->Velocity(idir,velx);
	} 	else if(iside==1) {
		molecule->Coordinate(idir,xmax[idir]-numeric_limits<REAL>::epsilon());
//-		molecule->Velocity(idir,-RND*velx);
		molecule->Velocity(idir,-velx);
	}  else {
		cerr<<"Inject: boundary not initialized\n";
		exit(1);
	}
	for(int i=0;i<DIM-1;i++) {
		int j=(idir+i+1)%DIM;
		REAL
			x=xmin[j]+RND*(xmax[j]-xmin[j]),
//-			v=(2.0*RND-1.0)*velx;
			v=vavx*Run::gauss();
		molecule->Coordinate(j,x);
		molecule->Velocity(j,v);
	}
	molecule->Temperature(temperature);
	GridContainer::put(molecule);
}
//! Domain constructor creates an instance of a domain  
Domain::Domain(char *filename) {
	//! get domain specs from XML configfile
	using namespace IO;
	using namespace Species;
	char *doctype=(char*)DOCTYPE;
	char *inpfilename=Run::configfile;
	//! read XML input file
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	for(int i=0;i<maxBoundaryTypes;i++)
		strcpy(boundaryName[i],"N/A");
	if(Run::option.verbose)printf("Domain::inpfilename=%s\n",inpfilename);
	xmlDocPtr doc;
	xmlNodePtr root;
	if (Run::option.verbose|Run::option.debug)
	{	printf
		(	"Reading %s\n",
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
	//! Loop over file records
	for 
	(	xmlNodePtr cur = root->xmlChildrenNode; 
		cur != NULL; cur = cur->next
	) {
		//! set xterm output mode:	
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"xterm"))) {
			char str[WORDLENGTH];
			xmlChar *key=xmlNodeListGetString(doc,cur->xmlChildrenNode, 1);
			strncpy(str,(char*)key,WORDLENGTH);
			xmlFree(key);
			int xterm=0;
			sscanf(str,"%d",&xterm);
			if(xterm!=0)Run::option.xterm=true;
			else Run::option.xterm=false;
			if(Run::option.verbose) {
				printf("Xterm output: %d\n",xterm);fflush(stdout);
			}
		}
		//! Read the number of molecules:	
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"molecules"))) {
			char str[WORDLENGTH];
			xmlChar *key=xmlNodeListGetString(doc,cur->xmlChildrenNode, 1);
			strncpy(str,(char*)key,WORDLENGTH);
			xmlFree(key);
			int nmolecules=0;
			sscanf(str,"%d",&nmolecules);
			if(Run::option.verbose) {
				printf("Maximum number of molecules: %d\n",nmolecules);fflush(stdout);
			}
			if(nmolecules>0) pool=new Pool<Molecule>(nmolecules);
			else pool=new Pool<Molecule>(1);
			molecules.init(nmolecules);
		}
		//! Read species:
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"species"))) {
			char title[MAXLINLEN];
			if(getCharAttr(cur,(char*)"title",title)==0) {
				fprintf
				(	stderr,"CAN'T GET %s FOR %s IN %s\nAborting\n",
					"name","domain",inpfilename
				); exit(1);
			}
			if(Run::option.verbose||Run::option.debug) printf("species: %s\n",title);
			xmlNodePtr tag=cur->xmlChildrenNode;
			//COUNT THE SPECIES:
			nspecies=0;
			while (tag != NULL) {	
				if ((!xmlStrcmp(tag->name, (const xmlChar *)"specie"))) 
					nspecies++;
				tag=tag->next;
			}
			if(Run::option.verbose) printf("Number of species: %d\n",nspecies);
			//ALLOCATE SPECIES:
			species=new Specie[nspecies+1];
			distribution=new int[nspecies];
			//FILL-IN THE SPECIES:
			tag=cur->xmlChildrenNode;
			for (int isp=0;tag != NULL&&isp<nspecies;tag = tag->next) {	
				Specie *specie=species+isp;
				if ((!xmlStrcmp(tag->name, (const xmlChar *)"specie"))) {
					xmlChar *key=xmlGetProp(tag,(const xmlChar *)"id");
					if(key!=NULL) {
						char id[WORDLENGTH];
						strncpy(id,(char*)key,WORDLENGTH);
						strncpy(specie->Id(),id,WORDLENGTH);
					}
					REAL mass=parseFloat(doc,tag,(char*)"mass");
					specie->Mass((REAL)mass);
					REAL size=parseFloat(doc,tag,(char*)"size");
					specie->Size((REAL)size);
					/*
						For monotomic molecules dof=3, which corresponds to 
						cp=(3/2+1)*R = 5/2*R
						since dof=2*cv/R=2*(cp/R-1)
						where R=BoltzmanConstant*AvogadroNumber 
						http://en.wikipedia.org/wiki/Specific_heat
						http://en.wikipedia.org/wiki/Heat_capacity#Heat_capacity
					*/
					REAL cp=parseFloat(doc,tag,(char*)"cp");// J/(mol*K) 
					if(cp<=2.5*GasConstant)cp=2.5*GasConstant;//can't be less than 3 dof
					specie->Cp(cp);
					if(Run::option.verbose){
						printf(
							"specie %d: name=%s, mass='%g', size='%g', cp=%g\n",
							isp,species[isp].Id(),species[isp].Mass(),species[isp].Size(),cp
						);fflush(stdout);
					}
					isp++;
				}
			}
			strcpy(species[VOIDSPECIE].Id(),"VOID\0");
			species[VOIDSPECIE].Size(TINY);//empty specie: to preclude poping up size in gui.cc
			species[VOIDSPECIE].Mass(0.0);//empty specie
			//allocate reactions:
			reactions=new Reaction[nspecies*nspecies];
			//! load reactions
			while(tag!=NULL) {
				if (!xmlStrcmp(tag->name, (const xmlChar *)"reaction")) {
					char reactants[WORDLENGTH],products[WORDLENGTH];
					xmlChar *key=xmlGetProp(tag,(const xmlChar *)"reactants");
					if(key==NULL){cerr<<"Keyword 'reactants' missing in reaction tag\n";exit(1);}
					strncpy(reactants,(char*)key,WORDLENGTH);
					key=xmlGetProp(tag,(const xmlChar *)"products");
					if(key==NULL){cerr<<"Keyword 'products' missing in reaction tag\n";exit(1);}
					strncpy(products,(char*)key,WORDLENGTH);
					REAL 
						activationEnergy=BoltzmannConstant/AtomicMassUnit*parseFloat(doc,tag,(char*)"activation"),
						probability=parseFloat(doc,tag,(char*)"probability"),
						// convert enthalpy from [kJ/mol] to [au*(nm/ns)^2]:
						enthalpy=1.0e3*parseFloat(doc,tag,(char*)"enthalpy")/(AvogadroNumber*AtomicMassUnit);
					//Parse reactants and products strings:
					char *r=reactants,*p=products,word[WORDLENGTH],*w=word;
					while(isspace(*r))r++; while(isspace(*p))p++;
					int ir0=0,ir1=0,ip0=0,ip1=0,j=0;
					for(w=word;w-word<WORDLENGTH-1&&!isspace(*r);r++,w++) *w=*r; *w='\0';
					for(ir0=0;ir0<nspecies;ir0++) if(strncmp(species[ir0].Id(),word,WORDLENGTH)==0) break;
					if(ir0==nspecies){cerr<<"First reactant not found for '"<<word<<"'\n";exit(1);}
					while(isspace(*r))r++;for(w=word;w-word<WORDLENGTH-1&&!isspace(*r);r++,w++) *w=*r; *w='\0';
					for(ir1=0;ir1<nspecies;ir1++) if(strncmp(species[ir1].Id(),word,WORDLENGTH)==0) break;
					if(ir1==nspecies){cerr<<"Second reactant not found for '"<<word<<"'\n";exit(1);}
					while(isspace(*p))p++;for(w=word;w-word<WORDLENGTH-1&&!isspace(*p);p++,w++) *w=*p; *w='\0';
					for(ip0=0;ip0<nspecies;ip0++) if(strncmp(species[ip0].Id(),word,WORDLENGTH)==0) break;
					if(ip0==nspecies){cerr<<"Product not found for '"<<word<<"'\n";exit(1);}
					while(*p!='\0'&&isspace(*p))p++;
					if(*p!='\0') {
						for(w=word;w-word<WORDLENGTH-1&&!isspace(*p);p++,w++) *w=*p; *w='\0';
						for(ip1=0;ip1<nspecies;ip1++) if(strncmp(species[ip1].Id(),word,WORDLENGTH)==0) break;
					} else ip1=VOIDSPECIE;
					int ireact0=ir0*nspecies+ir1, ireact1=ir1*nspecies+ir0;
					//Different reaction outcomes
					Reaction *react0=reactions+ireact0, *react1=reactions+ireact1;
					react0->Add(ip0,ip1,activationEnergy,probability,enthalpy);
					react1->Add(ip0,ip1,activationEnergy,probability,enthalpy);
					//Inverse reaction: include explicitly, since hard to handle different outcomes
				}
				tag=tag->next;
			}
			if(Run::option.verbose||Run::option.debug) {
				cout<<"Reactions between "<<nspecies<<" species:\n";
				for(int i=0;i<nspecies;i++)
				for(int j=0;j<nspecies;j++){
					Reaction *r=reactions+i*nspecies+j;
					Reaction::Outcome *o=r->First();
					if(o!=NULL)
					do{
						int ip0=o->Product(0), ip1=o->Product(1);
						REAL rate=o->Probability(), enthalpy=o->Enthalpy();
					if(rate>0.0) {
						cout<<i<<','<<j<<':'<<species[i].Id()<<" + "<<species[j].Id()<<" = "
							<<species[ip0].Id()<<" + "<<species[ip1].Id()
							<<"; rate="<<rate<<", ent="<<enthalpy<<endl;
							cout.flush();
					}//endif
					} while((o=r->Next())!=NULL);
				}
			}
		}//END species
		//! Read domain specs:
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))) {
			char name[MAXLINLEN];
			if(getCharAttr(cur,(char*)"name",name)==0) {
				fprintf
				(	stderr,"CAN'T GET name FOR domain IN %s\nAborting\n",
					inpfilename
				); exit(1);
			}
			if(Run::option.verbose) printf("Domain: %s\n",name);
			if(parseWord(doc,cur,(char*)"type",word)==0) {
				fprintf
				(	stderr,"Tag \"type\" is missing in %s\nAborting\n",
					inpfilename
				); exit(1);
			}
			if(Run::option.verbose){printf("Domain type: %s\n",word);fflush(stdout);}
			if(strcmp(word,"box")!=0) {
				fprintf(stderr,"ERROR: Only domain type 'box' is supported\n");
				exit(1);
			}
			//Only for box domain: list names as they
			//will appear in the xml file:
			strcpy(boundaryName[insideBoundary],"inside");
			strcpy(boundaryName[periodicBoundary],"periodic");
			strcpy(boundaryName[elasticBoundary],"elastic");
			strcpy(boundaryName[openBoundary],"open");
			for 
			(	xmlNodePtr next=cur->xmlChildrenNode;
				next != NULL;
				next=next->next
			) {
				if ((!xmlStrcmp(next->name, (const xmlChar *)"grid"))) {
					//Interaction acceleration scheme
					if(Run::option.verbose){printf("Particle partition grid:\n");fflush(stdout);}
					if(parseWord(doc,next,(char*)"cellsize",word)!=0) {
						REAL cellsize=0.0;
						sscanf(word,"%lg", &cellsize);
						GridContainer::cellsize=cellsize;
						if(Run::option.verbose) printf (	"\tiCell size: %g\n",cellsize);
					}
				} else//END grid 
				if ((!xmlStrcmp(next->name, (const xmlChar *)"energy"))) {//TO BE DEPRICATED
					if(Run::option.verbose){printf("Interaction Potential:\n");fflush(stdout);}
					if(parseWord(doc,next,(char*)"lengthscale",word)!=0) {
						REAL lengthscale=0.0,strength=0.0;
						sscanf(word,"%lg", &lengthscale);
						Potential::lengthscale(lengthscale);
						if(Run::option.verbose) printf (	"\tlengthscale: %g\n",lengthscale);
					}
					if(parseWord(doc,next,(char*)"strength",word)!=0) {
						REAL strength=0.0;
						sscanf(word,"%lg", &strength);
						Potential::strength(strength);
						if(Run::option.verbose) printf("\tstrength: %g\n",strength);
					}
				}//END energy 
				else//! Read domain bounds 
				if((!xmlStrcmp(next->name, (const xmlChar*)"bounds"))) {
					xmlChar *key=xmlNodeListGetString(doc,next->xmlChildrenNode,1);
					strncpy(word,(char*)key,MAXLINLEN);
					xmlFree(key);
					sscanf
					(	word,"%lg %lg %lg %lg %lg %lg",
						xmin,xmax,xmin+1,xmax+1,xmin+2,xmax+2
					);
					if(Run::option.verbose)
					printf
					(	"\tbounds: (%lg:%lg), (%lg:%lg), (%lg:%lg)\n",
						xmin[0],xmax[0],xmin[1],xmax[1],xmin[2],xmax[2]
					);
					//Check domain bounds
					bool adjusted=false;
					for(int i=0;i<DIM;i++) {
						if(xmin[i]>minBound(i)) {
							xmin[i]=minBound(i);
							adjusted=true;
						} else setMinBound(i,xmin[i]);
						if(xmax[i]<maxBound(i)) {
							xmax[i]=maxBound(i);
							adjusted=true;
						} else setMaxBound(i,xmax[i]);
					}
					if(adjusted) 
						printf
						(	"\tDomain bounds adjusted: (%lg:%lg), (%lg:%lg), (%lg:%lg)\n",
							xmin[0],xmax[0],xmin[1],xmax[1],xmin[2],xmax[2]
						);
					//set boundary areas:
					for(int idir=0;idir<DIM;idir++) {
						int i1=(idir+1)%DIM, i2=(idir+2)%DIM;
						REAL area=(xmax[i1]-xmin[i1])*(xmax[i2]-xmin[i2]);
						for(int iside=0;iside<2;iside++) {
							boundaries[idir][iside].area=area;
							if(Run::option.verbose)
								cout<<"Boundary ("<<idir<<','<<iside<<").area="<<area<<endl;
						}
					}
				}//END bounds
				else//! Read initial bulk properties
				if((!xmlStrcmp(next->name, (const xmlChar*)"bulk"))) {
					if(Run::option.verbose)printf("Domain bulk properties:\n");
					REAL temperature=300.0;//K: default bulk temperature
					if(parseFloat(doc,next,(char*)"temperature",temperature)!=0) {
						bulk.Temperature(temperature);
						if(Run::option.verbose)printf("\tbulk temperature=%g\n",temperature);
					}
					if(temperature<=0.0) {
						printf("Wrong domain temperature: %g\n",temperature);
						exit(1);
					}
					// Read bulk species:
					xmlNodePtr tag=next->xmlChildrenNode;
					while(tag!=NULL) {
						if(!xmlStrcmp(tag->name,(const xmlChar*)"specie")) {
							//Get Specie ID:
							char id[MAXLINLEN];
							if(getCharAttr(tag,(char*)"id",id)) {
								if(Run::option.verbose||Run::option.debug) {
									cout<<"\tBulk Specie id="<<id<<endl;cout.flush();
								}
								//Locate the element by ID:
								Specie *specie=NULL;
								for(int isp=0;isp<nspecies;isp++) {
									Specie *spc=species+isp;
									if(strcmp(id,spc->Id())==0){
										specie=spc;break;
									}
								}
								if(specie==NULL) {
									cerr<<"FAILED TO LOCATE SPECIE ID="<<id<<endl;
									exit(1);
								}
								//Append to the list of species
								Gas *gas=new Gas(specie);
								bulk.gases.append(gas);
								REAL density=0.0;
								if(parseFloat(doc,tag,(char*)"density",density)) 
									if(Run::option.verbose||Run::option.debug) {
										cout<<"\tDensity [kg/m^3]="<<density<<endl;
									}
								else 
									if(Run::option.verbose||Run::option.debug) {
										cout<<"Density not specified for specie "<<id<<", assumed 0\n";
									}
								gas->density=density;
							}	else {
								cerr<<"No Specie ID specified in "<<inpfilename<<endl;
								exit(1);
							}
						} 
						tag=tag->next;
					}
				}//END bulk
				else//! Read boundary specs
				if((!xmlStrcmp(next->name, (const xmlChar*)"boundary"))) {
					char id[MAXLINLEN];
					if(getCharAttr(next,(char*)"id",id)==0) {
						fprintf
						(	stderr,"CAN'T GET id FOR boundary IN %s\nAborting\n",
							inpfilename
						); exit(1);
					}
					int idir=-1,iside=-1;
					if(strcmp("left",id)==0) {idir=0;iside=0;}
					else if(strcmp("right",id)==0) {idir=0;iside=1;}
					else if(strcmp("bottom",id)==0) {idir=1;iside=0;}
					else if(strcmp("top",id)==0) {idir=1;iside=1;}
					else if(strcmp("front",id)==0) {idir=2;iside=0;}
					else if(strcmp("rear",id)==0) {idir=2;iside=1;}
					else {
						fprintf(stderr,"UNKNOWN BOUNDARY ID '%' in %s\n",id,inpfilename);
						exit(1);
					}
					Boundary *boundary=&boundaries[idir][iside];
					boundary->init(idir,iside,xmin,xmax,nspecies);
					if(parseWord(doc,next,(char*)"type",word)!=0) {
						int itype=0;
						p=word;
						while(isspace(*p))p++;
						for(itype=0;itype<maxBoundaryTypes;itype++) {
							char *name=boundaryName[itype];
							if(strncmp(p,name,strlen(name))==0) break;
						}
						if(itype==maxBoundaryTypes) {
							fprintf(stderr,"Can't assign boundary type '%s'\n",p);
							exit(1);
						}
						boundary->type=(BoundaryTypes)itype;
					}
					REAL temperature=-10.0;
					if(parseFloat(doc,next,(char*)"temperature",temperature)!=0) {
						boundary->temperature=temperature;
					}
					if(temperature<=0.0) boundary->adiabatic=true;
					else boundary->adiabatic=false;
					if(Run::option.verbose) {
						char *type=(char*)"isothermal";
						if(boundary->adiabatic) type=(char*)"adiabatic";
						printf(
							"Boundary %s (%d,%d): type=%d: %s, temperature=%g, heat-transfer: %s\n",
							id,idir,iside,(int)boundary->type,boundaryName[boundary->type],boundary->temperature,type
						);
					}
					xmlNodePtr tag=next->xmlChildrenNode;
					while(tag!=NULL) {
						if(!xmlStrcmp(tag->name,(const xmlChar*)"specie")) {
							//Get Specie ID:
							char id[MAXLINLEN];
							if(getCharAttr(tag,(char*)"id",id)) {
								if(Run::option.verbose||Run::option.debug) {
									cout<<"\tBoundary Specie id="<<id<<endl;cout.flush();
								}
								//Locate the element by ID:
								Specie *specie=NULL;
								for(int isp=0;isp<nspecies;isp++) {
									Specie *spc=species+isp;
									if(strcmp(id,spc->Id())==0){
										specie=spc;break;
									}
								}
								if(specie==NULL) {
									cerr<<"FAILED TO LOCATE SPECIE ID="<<id<<endl;
									exit(1);
								}
								//Append to the list of species
							//-	boundary->species.append(specie);
								Gas *gas=new Gas(specie);
								boundary->gases.append(gas);
								REAL density=0.0;
								if(parseFloat(doc,tag,(char*)"density",density)) 
									if(Run::option.verbose||Run::option.debug) {
										cout<<"\tDensity [kg/m^3]="<<density<<endl;
									}
								else 
									if(Run::option.verbose||Run::option.debug) {
										cout<<"Density not specified for specie "<<id<<", assumed 0\n";
									}
								gas->density=density;
							}	else {
								cerr<<"No Specie ID specified in "<<inpfilename<<endl;
								exit(1);
							}
						} 
						else if(!xmlStrcmp(tag->name,(const xmlChar*)"reaction")) {
							//! Read boundary reactions:
							//Get Reactant Specie ID:
							char reactant[WORDLENGTH],products[WORDLENGTH];
							xmlChar *key=xmlGetProp(tag,(const xmlChar *)"reactant");
							if(key==NULL){cerr<<"Keyword 'reactant' missing in boundary reaction tag\n";exit(1);}
							strncpy(reactant,(char*)key,WORDLENGTH);
							if(Run::option.verbose||Run::option.debug) {
									cout<<"\tBoundary Reactant ="<<reactant<<endl;cout.flush();
							}
							key=xmlGetProp(tag,(const xmlChar *)"products");
							if(key==NULL){cerr<<"Keyword 'products' missing in boundary reaction tag\n";exit(1);}
							strncpy(products,(char*)key,WORDLENGTH);
							if(Run::option.verbose||Run::option.debug) {
								cout<<"\tBoundary Products ="<<products<<endl;cout.flush();
							}
							REAL 
								activationEnergy=BoltzmannConstant/AtomicMassUnit*parseFloat(doc,tag,(char*)"activation"),
								probability=parseFloat(doc,tag,(char*)"probability"),/*!< probability for the 
								   formation of these products as opposed to different products
								   If only one pair of products exists the probability = 1
								*/
#ifdef ADSORPTION
								time=parseFloat(doc,tag,(char*)"time"),//!< time delay for the release from the surface
#endif
								enthalpy=parseFloat(doc,tag,(char*)"enthalpy");
								//Parse reactants and products strings:
							char *r=reactant,*p=products,word[WORDLENGTH],*w=word;
							while(isspace(*r))r++; while(isspace(*p))p++;
							int ir0=0,ip0=0,ip1=0,j=0;
							//Reactant:
							for(w=word;w-word<WORDLENGTH-1&&!isspace(*r);r++,w++) *w=*r; *w='\0';
							for(ir0=0;ir0<nspecies;ir0++) if(strncmp(species[ir0].Id(),word,WORDLENGTH)==0) break;
							if(ir0==nspecies){cerr<<"Boundary reactant not found for '"<<word<<"'\n";exit(1);}
							//Products:
							while(isspace(*p))p++;for(w=word;w-word<WORDLENGTH-1&&!isspace(*p);p++,w++) *w=*p; *w='\0';
							for(ip0=0;ip0<nspecies;ip0++) if(strncmp(species[ip0].Id(),word,WORDLENGTH)==0) break;
							if(ip0==nspecies){cerr<<"Boundary products not found for '"<<word<<"'\n";exit(1);}
							while(*p!='\0'&&isspace(*p))p++;
							if(*p!='\0') {
								for(w=word;w-word<WORDLENGTH-1&&!isspace(*p);p++,w++) *w=*p; *w='\0';
								for(ip1=0;ip1<nspecies;ip1++) if(strncmp(species[ip1].Id(),word,WORDLENGTH)==0) break;
							} else ip1=VOIDSPECIE;
#ifdef ADSORPTION
							boundary->reactions[ir0].Add(ip0,ip1,activationEnergy,probability,enthalpy,time);
#else
							boundary->reactions[ir0].Add(ip0,ip1,activationEnergy,probability,enthalpy);
#endif
						}
						tag=tag->next;
					}//end while
					if(Run::option.verbose||Run::option.debug) {
						cout<<"Reactions between "<<nspecies<<" boundary species:\n";
						for(int i=0;i<nspecies;i++){
							Reaction *r=boundary->reactions+i;
							Reaction::Outcome *o=r->First();
							if(o!=NULL) {
								do{
									int ip0=o->Product(0), ip1=o->Product(1);
									REAL rate=o->Probability(), enthalpy=o->Enthalpy();
									if(rate>0.0) {
										cout<<i<<':'<<species[i].Id()<<" = "
											<<species[ip0].Id()<<" + "<<species[ip1].Id()
											<<"; probability="<<rate<<", enthalpy="<<enthalpy<<endl;
											cout.flush();
									}
								} while((o=r->Next())!=NULL);
							}
						}
					}
				} 
			}//END next
		}//END if domain
	}//END cur
	xmlFreeDoc(doc);
	if(Run::option.verbose){printf("File %s parsed\n",inpfilename);fflush(stdout);}
	if(Run::option.restart)
		load(filename);// load molecules from the gzipped data file:
	else
		init();// initialize molecules according to bulk properties
#ifdef LOCAL
	GridContainer::init(
		xmin, xmax,
		//-1.5*Potential::Cutoff(), // grid-cell-size: heuristic
		GridContainer::cellsize,
		pool
	);
	if(Run::option.verbose||Run::option.debug) {
		if(Run::option.debug)cout<<"pool:"<<pool<<endl;cout.flush();
		cout<<"GridContainer::init:molecules:"<<molecules.number()<<endl;cout.flush();
	}
	if(molecules.number()>0) {
		molecules.goFirst();
		do {
			GridContainer::put(molecules.Current());
			molecules.goNext();
		}	while(!molecules.isFirst());
	}
#endif
}
void Domain::init() {
	bulk.init(xmin,xmax);
	double temperature=bulk.Temperature(),// get bulk temperature in K
		volume=(double)bulk.Volume()*1e-27;// volume of the bulk in m^3
	List<Gas> *gases=&bulk.gases;//! create local gas-array
	int imolecule,nmolecules=0;//total number of molecules and current molecule counter
	if(gases->number()>0) {
		gases->goFirst();
		do {//! Go over all gases
			Gas *gas=gases->Current();//! point to current gas
			Specie *specie=gas->specie;//! select current gas-specie
			/*! Estimate the average velocity and the number of molecules 
			 * of each specie to be injected into the bulk
			 */
			double // http://en.wikipedia.org/wiki/Kinetic_theory#Number_of_collisions_with_wall
				mu=(double)specie->Mass(),//! get molecule's mass in atomic units
				mass=(double)AtomicMassUnit*mu,//! convert molecule's mass to kg
//-				Vav=sqrt(8.0*BoltzmannConstant*temperature/(PI*mass)),//! estimate average velocity
				//See eq. (2)
				energy=BoltzmannConstant*temperature,//convert temperature to energy
				Vavx=sqrt(energy/mass),//! average velocity in one direction
				density=(double)gas->density,//! get density in kg/m^3
				totmass=density*volume;
			int	nmol=(int)round(totmass/mass);//! get the number of molecules in the bulk 
			//! Inject all molecules:
			if(Run::option.verbose) {
				printf("Domain::init:Injecting %d molecules of specie %s into the bulk\n",nmol,specie->Id());
			}
			for(int inj=0;inj<nmol;inj++) {
				Molecule *molecule = bulk.inject((int)(specie-Species::species),Vavx,&molecules);
				molecule->Temperature(temperature); // set internal energy
			}
			gases->goNext();
		}	while(!gases->isFirst());
	}
}
//! Initializing bulk: setting volume
void Bulk::init(REAL xmn[], REAL xmx[]) {
	volume=1.0;
	for(int i=0;i<DIM;i++) {
		xmin[i]=xmn[i];xmax[i]=xmx[i];
		volume*=xmax[i]-xmin[i];
	}
}
//! Injecting a specie into the bulk with a specified velocity (no registering in the grid)
Molecule *Bulk::inject(
	int ispecie, 
	REAL velx, // average velocity in one direction
	Collection<Molecule> *molecules
) {
	Molecule *molecule=molecules->append();
	molecule->Type(ispecie);
	//Select position: 
	for(int i=0;i<DIM;i++) {
		// random coordinates inside a box:
		REAL x=xmin[i]+RND*(xmax[i]-xmin[i]);
		molecule->Coordinate(i,x);
		//Select velocity from Gaussian: 
		REAL //- v=(2.0*RND-1.0)*velx;// uniform distribution
			v=velx*Run::gauss();// Gaussian distribution
		molecule->Velocity(i,v);
	}
	return molecule;
}
void Domain::load(char *filename) {
	uLong	nbuf=MAXLINLEN*sizeof(int);
	Byte *buf=(Byte*)calloc((uInt)nbuf,1);
	gzFile gzinp;
	if((gzinp = gzopen(filename, "rb"))==NULL)
	{	cerr << "CAN'T OPEN "<<filename<<endl;cerr.flush();
		exit(1);
	}
	if(Run::option.verbose||Run::option.debug)
		cout<<"Reading from "<<filename<<endl;
	while(true) {
		char	*p;
		gzgets(gzinp, (char*)buf,nbuf);
		if((char)buf[0]!='#')break;
		char *s=(char*)buf;
		if((p=strstr(s,"time"))!=NULL){
			p+=5;
			while(isspace(*p))p++;
			REAL starttime=0.0;
			sscanf(p,"%lg",&starttime);
			Run::time.current=starttime;
			if(Run::time.start>starttime) Run::time.start=starttime;
			if(Run::option.verbose) {cout<<"Start time="<<starttime<<endl;cout.flush();}
		}
		if((p=strstr(s,"output"))!=NULL){
			p+=6;
			while(isspace(*p))p++;
			int ioutput=0;
			sscanf(p,"%d",&ioutput);
			IO::ioutput=ioutput;//DDD +1;
			if(Run::option.verbose) {cout<<"Output dataset="<<IO::ioutput<<endl;cout.flush();}
		}
	}
	int imolecule,nmolecules=0;//total number of molecules and current molecule counter
	sscanf((char*)buf,"%d",&nmolecules);
	if(Run::option.verbose||Run::option.debug){cout<<"Number of molecules: "<<nmolecules<<endl;cout.flush();}
	REAL ymin[DIM],ymax[DIM];
	if(nmolecules>0) {
		for(int i=0;i<DIM;i++) {
			ymin[i]=LARGE;
			ymax[i]=-LARGE;
		}
	} else {
		for(int i=0;i<DIM;i++)ymin[i]=ymax[i]=0.0;
	}
	REAL maxsize=0.0;
//-		maxmass=0.0,
//-		minmass=LARGE;
	for(imolecule=0;imolecule<nmolecules&&!gzeof(gzinp);imolecule++) {
		char *s=(char*)buf;
		gzgets(gzinp, (char*)buf, nbuf);
		if(*s=='#') {imolecule--; continue;}
		Molecule *molecule=molecules.append();
		char id[WORDLENGTH];
		sscanf(s,"%s",&id);
		int type=-1;
		for(int i=0;i<nspecies;i++) {
			if(strncmp(id,species[i].Id(),WORDLENGTH)==0) {
				type=i;
				break;
			}
		}
		if(type<0) {
			cout<<"Unknown specie identifier "<<id<<endl;
			exit(1);
		}
		molecule->Type(type);
//-		REAL size=species[type].Size();
		REAL size=species[type].Size();
		if(maxsize<size)maxsize=size;
//-		REAL mass=species[type].Mass();
//-		if(maxmass<mass)maxmass=mass;
//-		if(minmass>mass)minmass=mass;
		while(!isspace(*++s));
		while(isspace(*++s));
		for (int i=0;i<DIM;i++) {//Read coordinates:
			REAL x;
			sscanf(s,"%lg",&x);
			molecule->Coordinate(i,x);
			if(ymin[i]>x)ymin[i]=x;
			if(ymax[i]<x)ymax[i]=x;
			while(!isspace(*++s));
			while(isspace(*++s));
		}
		for (int i=0;i<DIM;i++) {//Read velocity:
			REAL v;
			sscanf(s,"%lg",&v);
			molecule->Velocity(i,v);
			while(!isspace(*++s));
			while(isspace(*++s));
		}
		{	REAL e; //Read energy
			sscanf(s,"%lg",&e);
			molecule->InternalEnergy(e);
			while(!isspace(*++s));
			while(isspace(*++s));
		}
#ifdef DEBUG
		if(Run::option.debug) {
			cout<<imolecule<<": "<<id<<", x:";cout.flush();
			for(int i=0;i<DIM;i++) {
				cout<<" "<<molecule->Coordinate(i);
			}
			cout<<imolecule<<", v:";cout.flush();
			for(int i=0;i<DIM;i++) {
				cout<<" "<<molecule->Velocity(i);
			}
			cout<<endl;cout.flush();
		}
#endif
	}
	if(imolecule!=nmolecules) {
		cerr<<"Number of molecules specified ("<<nmolecules
		<<") does not agree with actually found ("<<imolecule<<")\n";
		cerr.flush();exit(1);
	}
	if(Potential::lengthscale()>0.0) Potential::Cutoff(Potential::lengthscale());
	else if(maxsize>0.0)Potential::Cutoff(2.0*maxsize);//DDD
	else { cerr<<"BAD LENGTHSCALE\n";cerr.flush();exit(1);}
	if(Run::option.debug) cout << molecules.number() 
		<< " molecules read. Interaction Cutoff = "<<Potential::Cutoff()<<endl;
	gzclose(gzinp);
}
int Domain::computeBounds(REAL x0[], REAL x1[]) {
	//computes domain bounds
	for(int i=0;i<DIM;i++) {
		x0[i]=LARGE;
		x1[i]=-LARGE;
	}
	molecules.goFirst();
	int nmolecules=0;
	do {
		Molecule *molecule=molecules.Current();
		REAL *x=molecule->Coordinates();
		for(int i=0;i<DIM;i++) {
			if(x0[i]>x[i])x0[i]=x[i];
			if(x1[i]<x[i])x1[i]=x[i];
		}
		molecules.goNext(); nmolecules++;
	}	while(!molecules.isFirst());
	return nmolecules;
}
//! The Domain solver: running <b>niter</b> iterations
int Domain::run(int niter) {
	static const REAL onethird=1.0/3.0;
	int *distribution=new int[nspecies];
	REAL currenttime=Run::time.current,
		starttime=Run::time.start,endtime=Run::time.end;
	if (Run::option.debug)
		cout << "Running "<< niter << " iterations ...";
	time_t	start_wtime=time((time_t) NULL);
	int iter=0;
	//! The main loop over niter iterations or till the end of time:
	for(;iter<niter&&currenttime<=endtime;iter++) {
		for(int i=0;i<nspecies;i++) {
			distribution[i]=0;
		}
		REAL 
			energy=0.0,// energy per DOF
			totalKineticEnergy=0.0,// kinetic energy per DOF
			totalInternalEnergy=0.0,// internal energy per DOF
			dtmin=LARGE; 

		//! Nested looping over all molecules:
		if(molecules.number()>0) {
			molecules.goFirst();
			do {
				Molecule	*molecule=molecules.Current();
				int type=molecule->Type();
				if(type!=VOIDSPECIE) {
					distribution[type]++;
#ifdef ADSORPTION
					if(molecule->ReleaseTime()<=currenttime) {
#endif
						molecule->Move();
						BoundaryTypes btype=boundary(molecule);
						type=molecule->Type();
						if(type==VOIDSPECIE) {
							molecules.remove();
							if(molecules.number()==0) break;
						} else {
							REAL
								cp=species[type].Cp(),
								dof=2.0*(cp*GasConstantInv-1.0),
								dofi=max(dof-3.0,0.0);//internal degrees of freedom
							REAL // energies per dof:
								internalEnergy=dofi>numeric_limits<REAL>::epsilon()?molecule->InternalEnergy()/dofi:(REAL)0.0,
								kineticEnergy=onethird*molecule->KineticEnergy();
							totalKineticEnergy+=kineticEnergy;
							totalInternalEnergy+=internalEnergy;
							energy+=0.5*internalEnergy+kineticEnergy;
#ifdef LOCAL	
							GridContainer::put(molecule);
#endif
							REAL 
								size=species[type].Size(),
								velocity=LENGTH(molecule->Velocity());
							if(velocity>SMALL)
							{	REAL dt=size/velocity;
								if(dt<dtmin) dtmin=dt;
							}
						}
#ifdef ADSORPTION
					}
#endif
				} 
				else {
					molecules.remove();
					if(molecules.number()==0) break;
				}
				molecules.goNext();
			}	while(!molecules.isFirst());
			Run::time.step=0.5*dtmin;// adjusting timestep 
			if(Run::option.xterm) {
				int ncols=6;//number of output colums on the xterm screen
				printf(
					"%c[HTime=%-9.3e Step=%-9.3e Memory=%d%% N=%d T=%-6.2fK E=%-9.3eJ K=%-9.3eJ I=%-9.3eJ%c[K\n",
					ESC,currenttime,Run::time.step,
					(int)(round(100.0*(REAL)molecules.number()/(REAL)molecules.maxnumber())),
					molecules.number(),
				//-	SQRTPIo8*AtomicMassUnit*energy/(BoltzmannConstant*(REAL)molecules.number()),ESC
				//-	AtomicMassUnit*kineticEnergy/(BoltzmannConstant*(REAL)molecules.number()),
				//-	AtomicMassUnit*internalEnergy/(BoltzmannConstant*(REAL)molecules.number()),
				//-	AtomicMassUnit*energy/(BoltzmannConstant*(REAL)molecules.number()),
					2.0*AtomicMassUnit*totalKineticEnergy/(BoltzmannConstant*(REAL)molecules.number()),
					AtomicMassUnit*energy,
					AtomicMassUnit*totalKineticEnergy,
					AtomicMassUnit*totalInternalEnergy,
					ESC
				);
				int ispecie=0,jspecie=0;
				while(ispecie<nspecies) {
					int icol=0;
					for(jspecie=ispecie;ispecie<nspecies;ispecie++){
						if(distribution[ispecie]>0) {
							if(icol++==ncols) break;
							printf("%%%-9s",species[ispecie].Id());
						}
					}
					if(ispecie<=nspecies)printf("%c[K\n",ESC);//ERASE END-OF-LINE
					icol=0;
					for(;jspecie<ispecie;jspecie++){
						if(distribution[jspecie]>0) {
							if(icol++==ncols) break;
							printf("%-9.5f ",100.0*(REAL)distribution[jspecie]/(REAL)molecules.number());
						}
					}
					if(ispecie<=nspecies)printf("%c[K\n",ESC);//ERASE END-OF-LINE
				}
				fflush(stdout);
			}
			interaction();
		} //endif molecules.number()>0
		injection();
		if(Run::time.output>0.0 && (int)floor(currenttime/Run::time.output)>IO::ioutput) {
			IO::ioutput=(int)floor(currenttime/Run::time.output);
			if(Run::option.verbose) cout<<"Writing output set "<<Run::outputname<<endl;
			save(Run::outputname);
		}
		currenttime+=Run::time.step;
		Run::time.current=currenttime;
	}
	clock_t end_wtime=time((time_t) NULL);
	REAL elapsed_time=difftime(end_wtime,start_wtime);
	if (Run::option.debug) 
		cout << " elapsed time: "<<elapsed_time<<"s, time per iteration: "<<elapsed_time/(REAL)niter<<"s\n";
	if(currenttime>endtime) currenttime=starttime+currenttime-endtime; // if looping
	if(currenttime<starttime) currenttime=endtime-starttime+currenttime; //return 0;
	Run::time.current=currenttime;
	return iter;
}
Domain::~Domain() {
	if(Run::option.debug) {cout<<"Deleting domain\n";cout.flush();}
//+	delete pool;
}
enum BoundaryTypes Domain::boundary(Molecule *a) {
	//PROCESS BOUNDARY CONDITIONS FOR PARTICLE a:
	REAL dx=0.0,
//-		radius=0.5*species[a->Type()].Size(),
		*x=a->Coordinates(),
		*v=a->Velocity();
	int idir=-1, iside=-1;//inside

	for(int i=0;i<DIM;i++) {
		if(x[i]<=xmin[i]) {idir=i;iside=0;break;}
		if(x[i]>=xmax[i]) {idir=i;iside=1;break;}
	}
	if(idir<0) return insideBoundary;
	dx=xmax[idir]-xmin[idir];
	Boundary *boundary=&boundaries[idir][iside];
	enum BoundaryTypes btype=boundary->type;
	switch(btype) {
		case periodicBoundary:
			// Periodic boundaries:
			if(iside==0) x[idir]+=dx;
			else         x[idir]-=dx;
//			for(int i=0;i<DIM;i++) {
//				REAL dx=xmax[i]-xmin[i];
//				if(x[i]<=xmin[i]) x[i]+=dx;
//				if(x[i]>=xmax[i]) x[i]-=dx;
//			}
			break;
		
		case elasticBoundary:// Elastic walls:
			v[idir]=-v[idir];
			if(iside==0) x[idir]=xmin[idir];
			else         x[idir]=xmax[idir];
//			{// Elastic wall collisions:
//				for(int i=0;i<DIM;i++) {
//					REAL dx=xmax[i]-xmin[i];
//					if(x[i]<=xmin[i]) {v[i]=-v[i];x[i]=xmin[i];}
//					else
//					if(x[i]>=xmax[i]) {v[i]=-v[i];x[i]=xmax[i];}
//				}
//			}
			break;
		case openBoundary://kill the molecule:
#ifdef LOCAL
			{
				Container<Molecule> *gridmolecules=GridContainer::nodes+a->GridCell();
				Ptr<Molecule> *ptr=a->getPtr();
				if(!gridmolecules->remove(ptr)) {//remove from the grid
					cout<<"Failed to remove node from the grid array\n";
					exit(1);
				}
			}
			a->putPtr(NULL);
#endif
			a->Type(VOIDSPECIE);// mark for removal
			return btype;
		default:
			cerr<<"Unknown BoundaryType "<<btype<<" encountered\n";
	}
	//Consider reactions:
	int ispecie=a->Type();
	Reaction *reaction=boundary->reactions+ispecie;
	Reaction::Outcome *outcome=reaction->First();
	if(outcome==NULL) return btype;
	REAL r=RND,//random number
		probability=0.0;
	bool reacted=false;
	do {
		probability+=outcome->Probability();
		if(r<probability) { //!< do reaction: A + B = 
			reacted=true;
			REAL 
				mAold=species[ispecie].Mass(),
				enthalpy=outcome->Enthalpy(),
#ifdef ADSORPTION
				time=outcome->Time(),//!< time for the reaction
#endif
				energy=a->KineticEnergy()+a->InternalEnergy()+enthalpy;
			int 
				typeA=outcome->Product(0),//!< new molecule A type
				typeB=outcome->Product(1);//!< new molecule B type
			if(typeA<0||typeA>=VOIDSPECIE) {
				cerr<<"Can't locate specie type "<<typeA<<" for boundary reaction\n";
				exit(1);
			}
			Specie *sA=species+typeA;
			REAL 
				mA=sA->Mass(),
				cpA=sA->Cp(),
				/*	Total degrees of freedom are distributed between 
 					the translational dof:
					doft=3
					and internal dof:
					dofi=dof-doft=2*cp/R-3 
					http://en.wikipedia.org/wiki/Heat_capacity#Heat_capacity
				*/
				dofA=2.0*(cpA/GasConstant-1.0);//=2*cv/R=degrees of freedom
			if(dofA<3.0)dofA=3.0;//check for physical consistency
			REAL	dofiA=dofA-3.0,
					vA=LENGTH(v),//velocity magnitude
					ieA=0.0;//internal energy of molecule A
			if(typeB!=VOIDSPECIE) {//create new molecule
				Specie *sB=species+typeB;
				REAL 
					mB=sB->Mass(),
					cpB=sB->Cp(),//heat capacity
					/*	Total degrees of freedom are distributed between 
 						the translational dof:
						doft=3
						and internal dof:
						dofi=dof-doft=2*cp/R-3 
						http://en.wikipedia.org/wiki/Heat_capacity#Heat_capacity
					*/
					dofB=2.0*(cpB/GasConstant-1.0);//=2*cv/R=degrees of freedom
				if(dofB<3.0)dofB=3.0;//check for physical consistency
				REAL	dofiB=dofB-3.0,//check for physical consistency
					ke=6.0/(dofA+dofB)*energy,//total kinetic energy
					ie=energy-ke,//total internal energy
					keA=0.5*ke,//kinetic energy of molecule A
					keB=keA,   //kinetic energy of molecule B
					dofi=dofiA+dofiB,
					ieB=0.0;//internal energy of B
				if(dofi>SMALL) {
					ieA=dofiA/dofi*ie;//internal energy of mol. A
					ieB=ie-ieA;//internal energy of molecule B
				}
				REAL
					m=mA+mB,
					vAnew=sqrt(2.0*keA/mA),
					vBnew=sqrt(2.0*keB/mB);
				if(vA>SMALL) {
					REAL c=sqrt(2.0*keA/mA)/vA;
					MULC(DIM,v,c);
				}
				Molecule *b=molecules.insert();
				//Rescale velocity to conserve energy:
				REAL 
					*u=b->Velocity(),
					*y=b->Coordinates();
				for(int i=0;i<DIM;i++){
					y[i]=x[i];
					u[i]=0.0;
				}
				REAL	d=sA->Size()+sB->Size();
				if(iside==0) {
					u[idir]=vBnew;
					y[idir]+=d+SMALL;//separate the two molecules 
				} else {
					u[idir]=-vBnew;
					y[idir]-=d+SMALL;//separate the two molecules 
				}
				b->InternalEnergy(ieB);
				b->Type(typeB);
#ifdef ADSORPTION
				if(time>SMALL) q->ReleaseTime(Run::time.current+time);// halt the molecule for time 'time'
#endif
			} else {//only one product
				//Rescale velocity to conserve energy:
				if(vA>SMALL) {
					REAL c=sqrt(2.0*energy/mA)/vA;
					MULC(DIM,v,c);
				}
			}
			a->InternalEnergy(ieA);
			a->Type(typeA);
#ifdef ADSORPTION
			if(time>SMALL) a->ReleaseTime(Run::time.current+time);// halt the molecule for time 'time'
#endif
		}
	}	while(reacted==false&&(outcome=reaction->Next())!=NULL);
	return btype;
}
//! Inject cross-boundary species:
void Domain::injection() {
	REAL dt=(REAL)Run::time.step*1e-9;//! get time-step in sec
	//! Go over all 6 boundaries:
	for(int idir=0;idir<DIM;idir++)
	for(int iside=0;iside<2;iside++) {
		Boundary *boundary=&boundaries[idir][iside];//! select a boundary
		REAL area=1.0e-18*(REAL)boundary->area,//! get boundry area in m^2
			temperature=boundary->temperature;//! get boundary temperature in K
		List<Gas> *gases=&boundary->gases;//! create local gas-array
		if(gases->number()>0) {
			gases->goFirst();
			do {//! Go over all gases
				Gas *gas=gases->Current();//! point to current gas
				Specie *specie=gas->specie;//! select current gas-specie
				/*! Use kinetic theory of gases to estimate collision frequency with the wall
 				 * and the number of molecules injected. This number can be fractional,
 				 * so split it into integer and fractional part and inject the fractional
 				 * part using probability.
 				 */
				REAL // http://en.wikipedia.org/wiki/Kinetic_theory#Number_of_collisions_with_wall
					mu=(REAL)specie->Mass(),//! get molecule's mass in atomic units
					mass=AtomicMassUnit*mu,//! convert molecule's mass to kg
					den=(REAL)gas->density,//! get density in kg/m^3
					energy=BoltzmannConstant*temperature,// energy per degree of freedom
					Vavx=sqrt(energy/mass),// average velocity in one direction
					velx=Vavx*fabs(Run::gauss()),// actual velocity
					freq=den*velx/(2.0*mass),//! estimate collision frequency per unit area
					totinj=freq*area*dt,//! calculate the number of collisions with the boundary during the past time step
					wholeinj=floor(totinj),//! get the whole integer number of molecules injected
					fracinj=totinj-wholeinj;//! calculate the fractional part of a molecule injected
				//! Inject all whole molecules:
				Molecule *molecule=NULL;
				int ninj=(int)wholeinj;
				for(int inj=0;inj<ninj;inj++) {
					boundary->Inject((int)(specie-Species::species),Vavx,velx,temperature,&molecules);
				}
				//! Probabilistically inject the fractional molecules:
				if(RND<=fracinj) {
					boundary->Inject((int)(specie-Species::species),Vavx,velx,temperature,&molecules);
				}
				gases->goNext();
			}	while(!gases->isFirst());
		}
	}
}
#ifdef LOCAL // ACCELERATION BY USING LOCAL INTERACTION MODEL
void Domain::interaction() {
	using namespace GridContainer;
	// Interaction of molecules a and b:
	int gridcell[DIM];
	REAL dx[DIM];
	for(int i=0;i<DIM;i++) dx[i]=xmax[i]-xmin[i];
	for(int ix=1;ix<ncells[0]-1;ix++) { gridcell[0]=ix;
	for(int iy=1;iy<ncells[1]-1;iy++) { gridcell[1]=iy;
	for(int iz=1;iz<ncells[2]-1;iz++) { gridcell[2]=iz;
		int icell=index(gridcell);
		Container<Molecule> *particles=nodes+icell;
		if(particles->number()==0) continue;
		particles->setFirstLast();
		do {
			Molecule *a=particles->First();
			if(a->Type()!=VOIDSPECIE) {
				Interaction status=missed;
				//Interact with the particles from the same GridContainer cell:
				if(particles->number()>1) {
					while (status==missed&&!particles->isLast()) {
						particles->goNext();
						Molecule *b=particles->Current();
						if(b->Type()!=VOIDSPECIE) {
							status=interact(a, b);
							if(status==annihalated) {
								if(b->Type()==VOIDSPECIE) {
									particles->remove(b->getPtr());
									b->putPtr(NULL);
									break;
								}
								else {
									cout<<"Molecule annihalated but not VOIDed!\n";cout.flush();exit(1);
								}
							}
						}
					}  // end while
				} //endif
				//Interact with the particles from adjacent GridContainer cells:
				int neibcell[DIM];
				for(int jx=-1;status==missed&&jx<2;jx++) { neibcell[0]=ix+jx;
				for(int jy=-1;status==missed&&jy<2;jy++) { neibcell[1]=iy+jy;
				for(int jz=-1;status==missed&&jz<2;jz++) { neibcell[2]=iz+jz;
					if(jx==0&&jy==0&&jz==0) continue;
					int jcell=index(neibcell);
					Container<Molecule> *neibs=nodes+jcell;
					if(neibs->number()==0) continue;
					neibs->goFirst();
					do {
						Molecule *b=neibs->Current();
						if(b->Type()!=VOIDSPECIE) {
							status=interact(a, b);
							if(status==annihalated) {
								if(b->Type()==VOIDSPECIE) {//Remove from the grid:
									neibs->remove(b->getPtr());
									b->putPtr(NULL);//the molecule is still listed
									//in the molecules collection, and will be
									//removed from there in the time iteration loop
									break;
								}
								else {
									cout<<"Molecule annihalated but not VOIDed!\n";cout.flush();exit(1);
								}
							}
						}
						neibs->goNext();
					}	while(status==missed&&!neibs->isFirst());
				}
				}
				}
			}//endif VOIDSPECIE
			particles->FirstNext();
		}	while(!particles->isFirstLast());
		particles->setFirst();
	}//end iz
	}//end iy
	}//end ix
}
#ifdef HARDBALLS
Interaction Domain::interact(Molecule *a, Molecule *b) {
	Interaction status=missed;
	int typeA=a->Type(), typeB=b->Type();
#ifdef DEBUG
	if(typeA==VOIDSPECIE||typeB==VOIDSPECIE){
		cerr<<"Wrong reactants species: \n"
		<<"\tA: type="<<typeA<<", id="<<species[typeA].Id()<<endl
		<<"\tB: type="<<typeB<<", id="<<species[typeB].Id()<<endl;
	}
#endif
	Specie 
		*spa=species+typeA,
		*spb=species+typeB;
	REAL
		*xa=a->Coordinates(), //coordinate of the point and it's periodic image
		*va=a->Velocity(),
		*xb=b->Coordinates();
	REAL 
		ma=spa->Mass(),
		ra=0.5*spa->Size(),//radius of a
		mb=spb->Mass(),
		rb=0.5*spb->Size(),//radius of b
		d[DIM],d0,d2=0.0;
	for(int i=0;i<DIM;i++)
	{	REAL r=xb[i]-xa[i];
		d[i]=r;
		d2+=r*r;
	}
	d0=sqrt(d2);
	if(d0<ra+rb) {//HARD BALL COLLISION: change velocities
		REAL	*vb=b->Velocity(),
		 	m=ma+mb,mi=1.0/m,//total mass and its inverse (for efficiency)
			d1[DIM],//unit vector in direction d
			w[DIM],//CM velocity
			ua[DIM],ub[DIM],//particle velocities
					// in CM coordinates
			una,unb;//lengths of ua and ub
		//Transfer to CM coordinates
		for(int k=0;k<DIM;k++)
		{	REAL 
				vak=(REAL)va[k],
				vbk=(REAL)vb[k],
				wk=mi*(ma*vak+mb*vbk);
			ua[k]=vak-wk;
			ub[k]=vbk-wk;
			w[k]=wk;
			d1[k]=d[k]/d0;
		}
		//Project CM velocities on distance vector:
		una=SCLP(ua,d1);
		unb=SCLP(ub,d1);
		if(una-unb>SMALL) {//to prevent molecule from locking on each-other
			REAL //Compute Energy:
				ua2old=SCLP(ua,ua),uaold=sqrt(ua2old),
				ub2old=SCLP(ub,ub),ubold=sqrt(ub2old),
				kineticEnergy=0.5*(ma*ua2old+mb*ub2old), //kinetic energy in CM
				internalEnergy=(REAL)a->InternalEnergy()+(REAL)b->InternalEnergy(),
				energy=kineticEnergy+internalEnergy;//total of two reacting molecules
			status=collided;
			//REVERSE VELOCITITES:
			for(int k=0;k<DIM;k++) {
				ua[k]-=2.0*d1[k]*una;
				ub[k]-=2.0*d1[k]*unb;
			}
			//Consider reaction
			Reaction *reaction=reactions+typeA*nspecies+typeB;
			Reaction::Outcome *outcome=reaction->First();
#ifdef DEBUG
			int typeAold=typeA, typeBold=typeB;
#endif
			if(outcome!=NULL) {// Reaction occurs
				REAL
					r=RND, //+Potential::invdist(RND),//add quantum uncertainty
					probability=0.0;
				do {
					probability+=outcome->Probability();
					if( r < probability &&
						 energy >= outcome->ActivationEnergy()
					) {//REACT: 
						typeA=outcome->Product(0);
						typeB=outcome->Product(1);
#ifdef DEBUG
						if(typeA==VOIDSPECIE) {//always erase b:DDD: Do we really need this?
							cout<<"typeA="<<typeA<<endl;///DDD
							typeA=typeB; typeB=VOIDSPECIE;
						}
#endif
						if(typeB!=VOIDSPECIE) {
							status=reacted;
						} else {//remove molecule b
							status=annihalated;
						}
						a->Type(typeA);
						spa=species+typeA;
						b->Type(typeB);
						spb=species+typeB;
						break;
					}
				}	while((outcome=reaction->Next())!=NULL);
			}
			REAL
				cpa=spa->Cp(),
				dofa=2.0*(cpa*GasConstantInv-1.0),//=2*cv/R=degrees of freedom
				cpb=spb->Cp(),
				dofb=2.0*(cpb*GasConstantInv-1.0),//=2*cv/R=degrees of freedom
				dof=dofa+dofb,//total degrees of freedom
				idof=0.0, //total internal degrees of freedom
				ridofa=0.0, ridofb=0.0;//percentile ratios of idofs for each molecule
			if(dof>6.0) { 
				idof=dof-6.0,//total internal degrees of freedom
				ridofa=dofa>3.0?(dofa-3.0)/idof:0.0;
				ridofb=1.0-ridofa;
			}
			if(status>=reacted) {//collision with changing mass:
				REAL 
					manew=spa->Mass(),
					mbnew=spb->Mass(),
					enthalpy=outcome->Enthalpy();
				energy+=enthalpy; // increase the energy by reaction enthalpy
#ifdef DEBUG
				if(fabs(ma+mb-manew-mbnew)>=SMALL) {
					cerr<<"Mass conservation violation in interaction:\n"
						<<species[typeAold].Id()<<'('<<species[typeAold].Mass()<<") + "
						<<species[typeBold].Id()<<'('<<species[typeBold].Mass()<<") = "
						<<ma+mb<<" != "
						<<spa->Id()<<'('<<spa->Mass()<<") + "
						<<spb->Id()<<'('<<spb->Mass()<<") = "
						<<manew+mbnew<<endl;
					exit(1);
				}
#endif
				if(status!=annihalated) {
					if(energy>0.0) {//Redistribute energy among DOFs:
						REAL
							//Translational:
							ke=6.0*energy/dof,//total kinetic energy of new molecules
							ie=energy-ke,//total internal energy of new molecules
							mratio=manew/mbnew,
							uanew=sqrt(2.0*(ke)/(manew*(1.0+mratio))),
							ubnew=uanew*mratio,
							uaratio=uanew/uaold,
							ubratio=ubnew/ubold;
						for(int i=0;i<DIM;i++) {
							ua[i]*=uaratio;
							ub[i]*=ubratio;
						}
						a->InternalEnergy(ridofa*ie);
						b->InternalEnergy(ridofb*ie);
					} else {//avoid negative energy: sqrt will fail
						REAL 
							ea=ridofa*enthalpy,
							eb=enthalpy-ea;
						a->InternalEnergy(a->InternalEnergy()+ea);
						b->InternalEnergy(b->InternalEnergy()+eb);
					}
				}	else {//fusion of two molecules with annihalation of b:
					for(int i=0;i<DIM;i++) { ua[i]=ub[i]=0.0; }
					a->InternalEnergy(energy);
					b->InternalEnergy(0.0);
				}
			}	else {//No reaction 
				REAL
					ke=6.0*energy/dof,//total kinetic energy of new molecules
					ie=energy-ke;//total internal energy of new molecules
				if(energy>0.0){//Redistribute internal energy among DOFs:
					REAL //Recalculate velocities to match ke:
						mratio=ma/mb,
						uanew=sqrt(2.0*(ke)/(ma*(1.0+mratio))),
						ubnew=uanew*mratio,
						uaratio=uanew/uaold,
						ubratio=ubnew/ubold;
					for(int i=0;i<DIM;i++) {
						ua[i]*=uaratio;
						ub[i]*=ubratio;
					}
					// Reset internal energies to match ie:
					a->InternalEnergy(ridofa*ie);
					b->InternalEnergy(ridofb*ie);
				} 
			} 
			//TRANSLATE TO LAB COORDINATES:
			for(int k=0;k<DIM;k++) {
				va[k]=(REAL)(ua[k]+w[k]);
				vb[k]=(REAL)(ub[k]+w[k]);
			}
		}//endif untilock
	}//endif d0<ra+rb
	return status;
}
#else //NOT-HARD-BALLS: Lennard-Jones:
bool Domain::interact(Molecule *a, Molecule *b) {
	bool collision=false;
	REAL 
//-		ra=0.5*a->Size(),//radius of a
		*xa=a->Coordinates(), //coordinate of the point and it's periodic image
//-		rb=0.5*b->Size(),//radius of ba
		*xb=b->Coordinates(),
		cutoffdistance=Potential::Cutoff(),
		d[DIM],d0,d2=0.0;
	for(int i=0;i<DIM;i++)
	{	REAL r=xb[i]-xa[i];
		d[i]=r;
		d2+=r*r;
	}
	d0=sqrt(d2);
	if(d0<cutoffdistance) {
		//Lennart-Jones interaction:
		// Compute force between the particles and change
		// particle velocity according to the force by
		// Newton's law: ua+=force/mass*dt
		REAL 
			dt=Run::time.step,
			f=Potential::force(d0),
			ma=species[a->Type()].Mass(),
			mb=species[b->Type()].Mass(),
			dva=f/ma*dt,
			dvb=f/mb*dt,
			*va=a->Velocity(),
			*vb=b->Velocity();
		for(int k=0;k<DIM;k++) {
			REAL r=d[k]/d0;
			va[k]+=dva*r; 
			vb[k]+=dvb*r;
		}
		collision=true;
	}
	return collision;
}
//ROLANDO: end of your function
#endif
#else // SIMPLE N^2 INTERACTION ALGORITHM:
void Domain::interaction() {
	// Interaction of molecules a and b:
	REAL dx[DIM];
	for(int i=0;i<DIM;i++) dx[i]=xmax[i]-xmin[i];
	molecules.setFirstLast();
	do {
		Molecule *a=molecules.First();
		REAL 
			ma=a->Mass(),
			ra=0.5*a->Size(),//radius of a
			*va=a->Velocity(),
			xa[2][DIM]; //coordinate of the point and it's periodic image
		for(int i=0;i<DIM;i++) {
			REAL x=a->Coordinate(i);
			xa[0][i]=x; // position of the point
			if(x-xmin[i]<0.5*dx[i])
				xa[1][i]=x+dx[i];//its periodic image
			else
				xa[1][i]=x-dx[i];
		}
		molecules.goFirstNext();
		while (!molecules.isLast()) {
			Molecule *b=molecules.Current();
			REAL 
				mb=b->Mass(),
				rb=0.5*b->Size(),//radius of ba
				*xb=b->Coordinates();
			for(int j=0;j<2;j++) {
				REAL d[DIM],d0,d2=0.0;
				for(int i=0;i<DIM;i++)
				{	REAL r=xb[i]-xa[j][i];
					d[i]=r;
					d2+=r*r;
				}
				d0=sqrt(d2);
				if(d0<ra+rb) { //HARD-BALLS COLLISION:
					REAL 	m=ma+mb,
						*vb=b->Velocity(),
						d1[DIM],//unit vector in direction d
						w[DIM],//CM velocity
						ua[DIM],ub[DIM],//particle velocities
								// in CM coordinates
						una,unb;//lengths of ua and ub
					//Transfer to CM coordinates
					for(int k=0;k<DIM;k++)
					{	w[k]=(ma*va[k]+mb*vb[k])/m;
						ua[k]=va[k]-w[k];
						ub[k]=vb[k]-w[k];
						d1[k]=d[k]/d0;
					}
					una=SCLP(ua,d1);
					unb=SCLP(ub,d1);
					if(una-unb>0.0) {
						for(int k=0;k<DIM;k++) {
							va[k]=ua[k]-2.0*d1[k]*una+w[k];
							vb[k]=ub[k]-2.0*d1[k]*unb+w[k];
						}
					}
				}
			} // end for j
			molecules.goNext();
		}  // end while
		molecules.FirstNext();
	}	while(!molecules.isFirstLast());
	molecules.setFirst();
}
#endif
void Domain::save(char *name) {
	char filename[MAXLINLEN];
	int nmolecules=molecules.number();
	if(nmolecules<=0) {
		if(Run::option.verbose) cout<< "No molecules saved\n";
		return;
	}
	sprintf(filename,"%s-%d.dat.gz",name,IO::ioutput);
	cout << "Saving task '"<<name<<"' into file: " << filename <<endl;
	gzFile out;
	if((out = gzopen(filename, "wb"))==NULL)
	{	cerr << "CAN'T OPEN "<<filename<<endl;cerr.flush();
		exit(1);
	}
	if(Run::option.verbose||Run::option.debug)
		cout<<"Writing "<<nmolecules<<" molecules to "<<filename<<" ... ";
	gzprintf(out,"# Data output %d at time %lg\n%d\n",IO::ioutput,Run::time.current,nmolecules);
	if(Run::option.verbose||Run::option.debug)cout<<"Number of molecules: "<<nmolecules<<endl;
	int imolecule=0;
	molecules.goFirst();
	do {
		Molecule *molecule=molecules.Current();
		if(molecule->Type()!=VOIDSPECIE) {
		gzprintf(out,"%s",species[molecule->Type()].Id());
		for(int i=0;i<DIM;i++) gzprintf(out,"\t%lg",molecule->Coordinate(i));
		for(int i=0;i<DIM;i++) gzprintf(out,"\t%lg",molecule->Velocity(i));
		gzprintf(out,"\t%lg\n",molecule->InternalEnergy());
#ifdef DEBUG
			if(Run::option.debug) {
				cout<<imolecule<<": size="<<species[molecule->Type()].Size()<<", x:";
				for(int i=0;i<DIM;i++) {
					cout<<" "<<molecule->Coordinate(i);
				}
				cout<<imolecule<<", v:";
				for(int i=0;i<DIM;i++) {
					cout<<" "<<molecule->Velocity(i);
				}
				cout<<endl;
			}
#endif
		}
		molecules.goNext();
		imolecule++;
	}	while(!molecules.isFirst());
	if(imolecule!=nmolecules) {
		cerr<<"Number of molecules specified ("<<nmolecules
		<<") does not agree with actually found ("<<imolecule<<")\n";
		cerr.flush();exit(1);
	}
	if(Run::option.verbose) cout << molecules.number() << " molecules written\n";
	gzclose(out);
//-	IO::ioutput++;
}

