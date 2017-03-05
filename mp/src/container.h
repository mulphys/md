template <class Element>
Container<Element>::Container() {
	nptrs=0;
	locked=false;
	current=first=last=NULL;
	for(int i=0;i<maxcounters;i++){
		current1[i]=first1[i]=NULL;
	}
}
template <class Element>
Container<Element>::~Container() {
	nptrs=0;
	current=first=last=NULL;
	for(int i=0;i<maxcounters;i++){
		current1[i]=first1[i]=NULL;
	}
}
template <class Element>
bool Container<Element>::checkPool(char *msg, Pool<Element> *pool) {
	if(nptrs>0) {
		current=first;
		do{
			if(pool->hook->next==current){cout<<"checkPool:hook->next=current="<<current<<endl;cout.flush();exit(1);}
			current=current->next;
		}while(current!=first);
	}
	return pool->check(msg);
}
template <class Element>
bool	Container<Element>::append0()
{	//expend the list by one ptr
	if(last==NULL) {
		last=GridContainer::pool->get();
		if(last==NULL) {
			cout<<"APPEND: CAN'T ALLOCATE FROM POOL\n";
			cout.flush();
			return false;
		}
		last->next=last->prev=last;
		last->element=NULL;
		current=first=last;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=last;
#endif
	} else {
		Ptr<Element> *lastnext=last->next;
		last->next=GridContainer::pool->get();
		if(last->next==NULL) {
			cout<<"APPEND: CAN'T ALLOCATE FROM POOL\n";cout.flush();
			return false;//failed to append: the memory is full
		}
		lastnext->prev=last->next;
		last->next->prev=last;
		last=last->next;
		last->next=lastnext;
		last->element=NULL;
		if(first!=lastnext){
			cout<<"ERROR IN append0: first="<<first<<" != last->next="<<lastnext<<endl;
			first=lastnext;
		}
	}
	nptrs++;
	return true;
}
template <class Element>
bool	Container<Element>::append()
{	//expend the list by one ptr
#ifdef OMP
	WAIT; locked=true;
#endif
	if(!append0())return false;
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return true;
}
template <class Element>
bool	Container<Element>::append(int icounter)
{	//expend the list by one ptr
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(first1[icounter]!=NULL)first=first1[icounter];
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	returned=append0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::append(Element *element)
{	//appends after last
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
#endif
	returned=append0(); 
	last->element=element;
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::append(Element *element,int icounter)
{	//appends after last
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(first1[icounter]!=NULL)first=first1[icounter];
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	returned=append0(); 
	last->element=element;
#ifdef OMP
	first1[icounter]=first;
	current1[icounter]=current;	
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
Ptr<Element>*	Container<Element>::append1(Element *element,int icounter)
{	//appends after last
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
//-	if(first1[icounter]!=NULL)first=first1[icounter];
//-	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	returned=append0(); 
	last->element=element;
#ifdef OMP
//-	first1[icounter]=first;
//-	current1[icounter]=current;	
//-	for(int i=0;i<maxcounters;i++) {
//-		if(first1[i]==NULL)first1[i]=first;
//-		if(current1[i]==NULL)current1[i]=current;
//-	}
	locked=false;
#endif
	return last;
}
template <class Element>
bool	Container<Element>::link0(Ptr<Element> *&ptr)
{	//links after last
	if(ptr==NULL) return false;
//-if(GridContainer::pool->hook->next==ptr){cout<<"link:ptr=hook->next="<<ptr<<", nptrs="<<nptrs<<endl;cout.flush();exit(1);}///DDD
	if(nptrs==0) {
		last=ptr;
		last->next=last->prev=last;
		current=first=last;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=last;
#endif
//-if(GridContainer::pool->hook->next==last){cout<<"Container::link:new:last=pool->hook->next="<<last<<endl;cout.flush();exit(1);}///DDD
	} else {
		//Link before last when nptrs>0:
		last->next->prev=ptr;
		ptr->next=last->next;
		ptr->prev=last;
		last->next=ptr;
		last=ptr;//make the new ptr last
	}
	nptrs++;
//-if(GridContainer::pool->hook->next==last){cout<<"Container::link:last=pool->hook->next="<<last<<endl;cout.flush();exit(1);}///DDD
//-#ifdef OMP
//-	for(int i=0;i<maxcounters;i++) {
//-		first1[i]=lastnext;
//-		current1[i]=current;
//-	}
//-#endif
	return true;
}
template <class Element>
bool	Container<Element>::link(Ptr<Element> *&ptr)
{	//links after last
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
#endif
	returned=link0(ptr);
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {///DDD
		if(first1[i]==NULL)first1[i]=first;///DDD
		if(current1[i]==NULL)current1[i]=current;///DDD	
	}///DDD
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::unlink0(Ptr<Element> *&ptr) {	// deletes current
	if(last==NULL||ptr==NULL||nptrs<=0) {
		cout<<"WARNING: Failed to unlink from the list of "<<nptrs<<": last="<<last<<", ptr="<<ptr<<endl;cout.flush();
//-if(first!=NULL)GridContainer::pool->put(first);
//-if(current!=NULL)GridContainer::pool->put(current);
//-if(ptr!=NULL)GridContainer::pool->put(ptr);
//-if(last!=NULL)GridContainer::pool->put(last);
//-first=current=last=NULL;nptrs=0;///DDD
//-for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=NULL;///DDD
//-return true;///DDD
		return false;
//-exit(1);///DDD
	}
	if(nptrs==1){
		if(ptr==first) {
			first=last=current=NULL;
#ifdef OMP
			for(int i=0;i<maxcounters;i++) first1[i]=current1[i]=NULL;
#endif
		} else {
			cerr<<"WARNING: CAN'T UNLINK POINTER: ptr="<<ptr<<", first="<<first<<", first->next="<<first->next<<", last="<<last<<", last->prev="<<last->prev<<", nptrs="<<nptrs<<endl;cerr.flush();
//-if(first!=NULL)GridContainer::pool->put(first);
//-if(current!=NULL)GridContainer::pool->put(current);
//-if(ptr!=NULL)GridContainer::pool->put(ptr);
//-if(last!=NULL)GridContainer::pool->put(last);
//-first=current=last=NULL;nptrs=0;///DDD
//-for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=NULL;///DDD
//-return true;///DDD
//+		last=ptr=NULL;nptrs=0;///DDD
//-return true;///DDD
			return false;//or dump data and abort
//-exit(1);///DDD
		}
	} else {
		if(ptr==current) current=current->next;
		if(ptr==first)first=first->next;
#ifdef OMP
		for(int i=0;i<maxcounters;i++){
			if(ptr==current1[i])current1[i]=current1[i]->next;
			if(ptr==first1[i])first1[i]=first1[i]->next;
		}
#endif
		if(ptr==last)last=last->prev;
		ptr->next->prev=ptr->prev;
		ptr->prev->next=ptr->next;
	}
	nptrs--;
	return true;
}
template <class Element>
bool	Container<Element>::unlink(Ptr<Element> *&ptr) {	// deletes current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
#endif
	returned=unlink0(ptr);
#ifdef OMP
for(int i=0;i<maxcounters;i++){///DDD
if(current1[i]==NULL)current1[i]=current;///DDD
if(first1[i]==NULL)first1[i]=first;}///DDD
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::insert0()
{
	if(first==NULL) {	//inserts first 
		first=GridContainer::pool->get(); //-new Ptr<Element>;
		if(first==NULL) {
			cout<<"FAILED TO INSERT INTO CONTAINER\n";cout.flush();
			return false;//failed to insert: the memory is full
		}
		first->next=first->prev=first;
		first->element=NULL;
		current=last=first;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=last;
#endif
	} else {
		//Insert at current:
		Ptr<Element> *newptr=GridContainer::pool->get();
		if(newptr==NULL) {
			cout<<"FAILED TO INSERT INTO CONTAINER\n";cout.flush();
			return false;//failed to insert: the memory is full
		}
		current->next->prev=newptr;
		newptr->next=current->next;
		newptr->prev=current;
		current->next=newptr;
		current=current->next;
		current->element=NULL;
		last=current;
	}
	nptrs++;
	return true;
}
template <class Element>
bool	Container<Element>::insert()
{	//inserts after current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
#endif
	returned=insert0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::insert(int icounter)
{	//inserts after current1[icounter]
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(first1[icounter]!=NULL)first=first1[icounter];
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	returned=insert0();
#ifdef OMP
	first1[icounter]=first;
	current1[icounter]=current;	
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::insert(Element *element)
{	//inserts after current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
#endif
	returned=insert0();
	current->element=element;
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::insert(Element *element,int icounter)
{	//inserts after current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(first1[icounter]!=NULL)first=first1[icounter];
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	returned=insert0();
	current->element=element;
#ifdef OMP
	first1[icounter]=first;
	current1[icounter]=current;	
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
Ptr<Element>*	Container<Element>::insert1(Element *element,int icounter)
{	//inserts after current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(first1[icounter]!=NULL)first=first1[icounter];
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	returned=insert0();
	current->element=element;
#ifdef OMP
	first1[icounter]=first;
	current1[icounter]=current;	
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return current;
}
template <class Element>
bool	Container<Element>::remove0() {	// deletes current, but does not delete the element
	if(current==NULL||nptrs<=0) {
		cout<<"CAN'T REMOVE FROM AN EMPTY CONTAINER OF "<<nptrs<<" ELEMENTS\n";cout.flush();
		return false;//failed to remove: the list is empty
	}
	if(nptrs==1) {
		GridContainer::pool->put(current);
		current=first=last=NULL;
#ifdef OMP
		for(int i=0;i<maxcounters;i++) first1[i]=current1[i]=NULL;
#endif
	} else {
		if(current==last)last=last->prev;
		if(current==first)first=first->next;
		Ptr<Element> *old=current;
		current=current->prev;
		GridContainer::pool->put(old);
	}
	nptrs--;
	return true;
}
template <class Element>
bool	Container<Element>::remove() {// deletes current, but does not delete the element
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	for(int i=0;i<maxcounters;i++){
		if(first1[i]==first)first1[i]=NULL;
		if(current1[i]==current)current1[i]=NULL;
	}
#endif
	returned=remove0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) {
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::remove(int icounter) {	// deletes current1[icounter], but does not delete the element
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(first1[icounter]!=NULL)first=first1[icounter];
	if(current1[icounter]!=NULL)current=current1[icounter];
	for(int i=0;i<maxcounters;i++){
		if(first1[i]==first)first1[i]=NULL;
		if(current1[i]==current)current1[i]=NULL;
	}
#endif
	returned=remove0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++){
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::remove0(Ptr<Element> *&ptr) {// deletes current
	if(current==NULL||ptr==NULL||nptrs<=0) {
		cout<<"Failed to remove from the list: the list is empty\n";cout.flush();
		return false;
	}
//-if(GridContainer::pool->hook->next==ptr){cout<<"Container::remove("<<ptr<<"):hook->next=ptr="<<ptr<<endl;cout.flush();exit(1);}///DDD
	if(nptrs==1){
		if(ptr==first) {
			GridContainer::pool->put(ptr);
			first=last=current=NULL;
#ifdef OMP
			for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=NULL;
#endif
		} else {
			cout<<"CAN'T REMOVE ptr="<<ptr<<" FROM THE CONTAINER\n";cout.flush();
			return false;//or dump and abort
		}
	} else {
		if(ptr==last)last=last->prev;
		if(ptr==current) current=current->next;
		if(ptr==first)first=first->next;
#ifdef OMP
		for(int i=0;i<maxcounters;i++) {
			if(ptr==current1[i])current1[i]=current1[i]->next;
			if(ptr==first1[i])first1[i]=first1[i]->next;
		}
#endif
		GridContainer::pool->put(ptr);
	}
//-	if(nptrs==0) {
//-		cout<<"INTERNAL ERROR: Negative pointer: nptrs=%d\n";cout.flush();
//-		exit(1);
//-	}
	nptrs--;
//+	current=ptr; //fast but unsafe
//+	remove();
	return true;
}
template <class Element>
bool	Container<Element>::remove(Ptr<Element> *ptr) {	// deletes current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
#endif
	returned=remove0(ptr);
#ifdef OMP
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Container<Element>::remove(Ptr<Element> *ptr,int icounter) {	// deletes current
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(icounter>=0&&icounter<maxcounters){
		if(current1[icounter]!=NULL)current=current1[icounter];
		if(first1[icounter]!=NULL)first=first1[icounter];
	}
	for(int i=0;i<maxcounters;i++){
		if(first1[i]==first)first1[i]=NULL;
		if(current1[i]==current)current1[i]=NULL;
	}
#endif
	returned=remove0(ptr);
#ifdef OMP
	for(int i=0;i<maxcounters;i++){
		if(first1[i]==NULL)first1[i]=first;
		if(current1[i]==NULL)current1[i]=current;
	}
	locked=false;
#endif
	return returned;
}
