#include <stdlib.h>
#include <cstdio>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "def.h"
#include "run.h"
#include "io.h"

namespace IO {
	int	ioutput=0,//output counter 
		nxwdump=0;//Window dump counter
	void	initxwd(int i)
	{
		nxwdump=i;
	}
	void	xwd(char *windowname)
	{
		int	sysreturned;
		char
			*s,
			dumpfilename[MAXLINLEN],
			command[MAXLINLEN];
		for (s=windowname; !isalpha(*s); s++);
		sprintf(dumpfilename,"%s%03d.xwd",s,nxwdump++);
		sprintf(command,"xwd -name %s -out %s",windowname,dumpfilename);
		if (Run::option.verbose||Run::option.debug)
		{	printf("Dumping window to '%s' ... ",dumpfilename);fflush(stdout);
		}
		sysreturned=system(command);
		if (Run::option.verbose||Run::option.debug) 
			printf("done\n");
		if (sysreturned!=0)
		{	fprintf(stderr,"Window dump to %s failed with status %d\n",sysreturned);
			fflush(stderr);
		}
	}
int getCharAttr(xmlNodePtr cur, char *attr, char *result) {
	xmlChar *key;
	key = xmlGetProp(cur, (const xmlChar *)attr);
	if(key!=NULL)
	{	strncpy(result,(char*)key,MAXLINLEN);
		xmlFree(key);
		return 1;
	}
	if(Run::option.verbose||Run::option.debug)
	{	fprintf
		(	stderr,
			"\tAttribute %s not found\n",
			attr
		);
	}
	return 0;
}
int getIntAttr(xmlNodePtr cur, char *attr) 
{
	xmlChar *key;
	key = xmlGetProp(cur, (const xmlChar *)attr);
	if(key!=NULL)
	{	int result;
		result=atoi((char*)key);
		xmlFree(key);
		return result;
	}
	if(Run::option.verbose||Run::option.debug)
	{	fprintf
		(	stderr,
			"\tAttribute %s not found\n",
			attr
		);
	}
	return -1;
}
int getIntAttr
(	xmlDocPtr doc, 
	char *tag, 
	char *keyname, 
	char *key, 
	char *attr
) 
{
	xmlNodePtr cur;
///	xmlChar *key;
	cur = xmlDocGetRootElement(doc);
	cur = cur->xmlChildrenNode;
	if(Run::option.verbose) 
	{	printf
		(	"getIntAttr:tag=%s, keyname=%s, key=%s, attr=%s\n",
			tag,keyname,key,attr
		);	fflush(stdout);
	}
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)tag)))
		{	char keyval[MAXLINLEN];
			if(getCharAttr(cur,keyname,keyval)!=0)
			{	if(strcmp(keyval,key)==0)
				{	int	value=0;
					///val=getIntAttr(cur,attr);
					value=parseInt(doc,cur,attr);
					return value;
				}
			}
		}
		cur=cur->next;
	}
	if(Run::option.verbose||Run::option.debug) {
		fprintf
		(	stderr,
			"\tTag '%s' not found\n",
			keyname
		); fflush(stderr);
	}
	return -1;
}
int parseWord(xmlDocPtr doc, xmlNodePtr cur, char *keyword, char *result) {
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{	key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			strncpy(result,(char*)key,MAXLINLEN);
			xmlFree(key);
			return 1;
		}
		cur = cur->next;
	}
	if(Run::option.verbose||Run::option.debug) {
		fprintf
		(	stderr,
			"\tTag '%s' not found\n",
			keyword
		); fflush(stderr);
	}
	return 0;
}
int parseInt(xmlDocPtr doc, xmlNodePtr cur, char *keyword) 
{
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{	xmlChar *key;
			int keyval;
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			keyval=(int)atoi((char*)key);
			xmlFree(key);
			return keyval;
		}
		cur = cur->next;
	}
	if(Run::option.verbose||Run::option.debug) {
		fprintf(stderr,"Integer data for tag '%s' not found\n",keyword);
		fflush(stderr);
	}
	return -1;
}
int parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword, REAL &result) {
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			//printf("%s: %s\n", keyword,key);///DDD
			if(key!=NULL) result=(REAL)atof((char*)key);
			xmlFree(key);
			return 1;
		}
		cur = cur->next;
	}
	if(Run::option.verbose||Run::option.debug) {
		fprintf(stderr,"Float data for tag '%s' not found\n",keyword);
		fflush(stderr);
	}
	return 0;
}
double parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword) {
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword))) 
		{	double result;
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			//printf("%s: %s\n", keyword,key);///DDD
			result=(double)atof((char*)key);
			xmlFree(key);
			return result;
		}
		cur = cur->next;
	}
	if(Run::option.verbose||Run::option.debug) {
		fprintf(stderr,"Float data for tag '%s' not found\n",keyword);
		fflush(stderr);
	}
	return 0.0;
}
void	getTimeXML(char inpfilename[])
{//?	extern double parseFloat(xmlDocPtr doc, xmlNodePtr cur, char *keyword);
	using namespace Run;
	int	level=0,timefound=0;
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	const char *doctype=DOCTYPE;
	xmlDocPtr doc;
	xmlNodePtr cur;
	//Default times
	Run::time.start=0.0;
	Run::time.end=1.0;
	Run::time.step=Run::time.step0=0.01;
	doc = xmlParseFile(inpfilename);
	if (doc == NULL ) 
	{	fprintf(stderr,"XML parser failed in %s\n",inpfilename);
		return;
	}
	cur = xmlDocGetRootElement(doc);
	if (cur == NULL) 
	{	fprintf(stderr,"Empty document: %s\n",inpfilename);
		xmlFreeDoc(doc);
		return;
	}
	if (xmlStrcmp(cur->name, (const xmlChar *) doctype)) 
	{	fprintf(stderr,"document of the wrong type, root node != %s in %s\n",doctype,inpfilename);
		xmlFreeDoc(doc);
		return;
	}
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)"time")))
		{	Run::time.start=parseFloat(doc,cur,(char*)"start");
			Run::time.end=parseFloat(doc,cur,(char*)"end");
			Run::time.step=parseFloat(doc,cur,(char*)"step");
			Run::time.output=parseFloat(doc,cur,(char*)"output");
		}
		cur = cur->next;
	}
	xmlFreeDoc(doc);
}
void	getTime(char inpfilename[])
{	using namespace Run;
	if(strstr(inpfilename,".xml\0")!=NULL)
		getTimeXML(inpfilename);
	else {
		fprintf(stderr,"File %s should have extention 'xml'\n"); exit(1);
	}
	if (Run::option.verbose)
		printf("Time: start=%g, end=%g, step=%g, output=%g\n",Run::time.start,Run::time.end,Run::time.step,Run::time.output);
}
int	getIter(char inpfilename[])
{	using namespace Run;
	if(strstr(inpfilename,".xml\0")==NULL) {
		fprintf(stderr,"File %s should have extention 'xml'\n"); exit(1);
	}
	int	niter=0,level=0;
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	const char *doctype=DOCTYPE;
	xmlDocPtr doc;
	xmlNodePtr cur;
	doc = xmlParseFile(inpfilename);
	if (doc == NULL ) 
	{	fprintf(stderr,"XML parser failed in %s\n",inpfilename);
		return 0;
	}
	cur = xmlDocGetRootElement(doc);
	if (cur == NULL) 
	{	fprintf(stderr,"Empty document: %s\n",inpfilename);
		xmlFreeDoc(doc);
		return 0;
	}
	if (xmlStrcmp(cur->name, (const xmlChar *) doctype)) 
	{	fprintf(stderr,"document of the wrong type, root node != %s in %s\n",doctype,inpfilename);
		xmlFreeDoc(doc);
		return 0;
	}
	cur = cur->xmlChildrenNode;
	while (cur != NULL) 
	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)"iterations"))) {
			xmlChar *key;
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			niter=(int)atoi((char*)key);
			xmlFree(key);
			break;
		}
		cur = cur->next;
	}
	xmlFreeDoc(doc);
	if (Run::option.verbose)
		printf("Read number of iterations %d from file %s\n",niter,inpfilename);
	return niter;
}
//-int	countObjectsXML
//-(	char	*keyword,
//-	char	*filename
//-)
//-{	int	count;
//-	xmlDocPtr doc;
//-	xmlNodePtr cur;
//-	doc = xmlParseFile(filename);
//-	if (doc == NULL ) 
//-	{	fprintf(stderr,"Document not parsed successfully in %s\n",filename);
//-		return 0;
//-	}
//-	cur = xmlDocGetRootElement(doc);
//-	if (cur == NULL) 
//-	{	fprintf(stderr,"empty document in %s\n",filename);
//-		xmlFreeDoc(doc);
//-		return 0;
//-	}
//-	if (xmlStrcmp(cur->name, (const xmlChar *) DOCTYPE)) 
//-	{	fprintf
//-		(	stderr,"document of the wrong type, root node != %s in %s\n",
//-			DOCTYPE,filename
//-		);
//-		xmlFreeDoc(doc);
//-		return 0;
//-	}
//-	cur = cur->xmlChildrenNode;
//-	count=0;
//-	while (cur != NULL) 
//-	{	if ((!xmlStrcmp(cur->name, (const xmlChar *)keyword)))
//-			count++;
//-		cur = cur->next;
//-	}
//-	xmlFreeDoc(doc);
//-	return count;
//-}
//-int	countObjectsCFG
//-(	char	*keyword,
//-	char	*filename
//-)
//-{	int	count,wordlen=strlen(keyword);
//-	char	*s,line[MAXLINLEN];
//-	FILE	*inp;
//-	OPENREAD(filename,inp);
//-	count=0;
//-	while ((s=fgets(line,MAXLINLEN,inp))!=NULL)
//-	{
//-		while(isspace(*s))s++;
//-		if (memcmp(s,keyword,wordlen)==0&&isspace(*(s+wordlen))) count++;
//-	}
//-	fclose(inp);
//-	return count;
//-}
//-int	countObjects
//-(	char	*keyword,
//-	char	*filename
//-)
//-{	int count;
//-	if(strstr(filename,".xml\0")!=NULL)
//-		count=countObjectsXML(keyword,filename);
//-	else
//-		count=countObjectsCFG(keyword,filename);
//-	if (Run::option.verbose)
//-		printf("Number of domains: count=%d\n",count);
//-	return count;
//-}
} //END namespace Input
