template <class Element>
Container<Element>::Container() {
	nptrs=0;
	current=first=last=NULL;
}
template <class Element>
Container<Element>::~Container() {
	nptrs=0;
	current=first=last=NULL;
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
bool	Container<Element>::append()
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
		nptrs=1;
//-if(last==GridContainer::pool->hook->next){cout<<"Container::append:last=hook->next="<<last<<endl;;cout.flush();exit(1);}///DDD
		return true;
	}
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
	nptrs++;
	first=lastnext;
//-if(first->next==NULL){cout<<"Container::append:GOTCHA!!!\n";cout.flush();exit(1);}///DDD
	return true;
}
template <class Element>
bool	Container<Element>::append(Element *element)
{	//appends after last
	if(append()==false) return false; 
	last->element=element;
	return true;
}
template <class Element>
bool	Container<Element>::link(Ptr<Element> *ptr)
{	//links after last
	if(ptr==NULL) return false;
//-if(GridContainer::pool->hook->next==ptr){cout<<"link:ptr=hook->next="<<ptr<<", nptrs="<<nptrs<<endl;cout.flush();exit(1);}///DDD
	if(nptrs==0) {
		last=ptr;
		last->next=last->prev=last;
		current=first=last;
		nptrs=1;
//-if(GridContainer::pool->hook->next==last){cout<<"Container::link:new:last=pool->hook->next="<<last<<endl;cout.flush();exit(1);}///DDD
		return true;
	}
	//Link before last when nptrs>0:
	last->next->prev=ptr;
	ptr->next=last->next;
	ptr->prev=last;
	last->next=ptr;
	last=ptr;//make the new ptr last
	nptrs++;
//-if(GridContainer::pool->hook->next==last){cout<<"Container::link:last=pool->hook->next="<<last<<endl;cout.flush();exit(1);}///DDD
	return true;
}
template <class Element>
bool	Container<Element>::unlink(Ptr<Element> *ptr) {	// deletes current
	if(last==NULL||ptr==NULL) {
		cout<<"Failed to unlink from the list";
		if(last==NULL) cout<<": the list is empty";
		cout<<endl;cout.flush();
		return false;
	}
//-cout<<"unlink:nptrs="<<nptrs<<", ptr="<<ptr<<", pool="<<GridContainer::pool<<", hook="<<GridContainer::pool->hook;cout.flush();///DDD
//-cout<<", hook->next="<<GridContainer::pool->hook->next<<endl;cout.flush();///DDD
//-if(GridContainer::pool->hook->next==ptr){cout<<"unlink:ptr=hook->next="<<ptr<<", nptrs="<<nptrs<<endl;cout.flush();exit(1);}///DDD
	if(nptrs==1){
//-cout<<"unlink nptrs=1: first="<<first<<", last="<<last<<", current="<<current;cout.flush();///DDD
//-cout<<"first->next="<<first->next<<endl;cout.flush();///DDD
		if(ptr==first) {
			first=last=current=NULL;
			nptrs=0;
			return true;
		} else {
			cout<<"CAN'T UNLINK POINTER: ptr="<<ptr<<", first="<<first<<", first->next="<<first->next<<", last="<<last<<", nptrs="<<nptrs<<endl;cout.flush();
		return false;//or dump data and abort
		}
	}
	if(ptr==current) current=current->next;
	if(ptr==last)last=last->prev;
	if(ptr==first)first=first->next;
	ptr->next->prev=ptr->prev;
	ptr->prev->next=ptr->next;
	nptrs--;
	return true;
}
template <class Element>
bool	Container<Element>::insert()
{	//inserts after current
	if(first==NULL) {
		first=GridContainer::pool->get(); //-new Ptr<Element>;
		if(first==NULL) {
			cout<<"FAILED TO INSERT INTO CONTAINER\n";cout.flush();
			return false;//failed to insert: the memory is full
		}
		first->next=first->prev=first;
//-cout<<"\nContainer:"<<nptrs<<": insert first->next="<<first->next<<endl;cout.flush();///DDD
		first->element=NULL;
		current=last=first;
		nptrs=1;
//-if(GridContainer::pool->hook->next==first){cout<<"Container::insert:hook->next==first="<<first<<endl;cout.flush();exit(1);}///DDD
		return true;
	}
	//Insert in front of the current:
	Ptr<Element> *newptr=GridContainer::pool->get();
//-cout<<"\nContainer:"<<nptrs<<": insert current->next="<<newptr<<endl;cout.flush();///DDD
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
	nptrs++;
//-if(newptr==GridContainer::pool->hook->next){cout<<"Container::insert:hook->next=newptr="<<newptr<<endl;cout.flush();exit(1);}///DDD
//-if(first->next==NULL){cout<<"Container::insert:GOTCHA!!!!!\n";cout.flush();exit(1);}///DDD
	return true;
}
template <class Element>
bool	Container<Element>::insert(Element *element)
{	//inserts after current
//-cout<<"Container::insert:element="<<element;cout.flush();///DDD
	if(!insert())return false;
	current->element=element;
//-cout<<", current="<<current;cout.flush();///DDD
	return true;
}
template <class Element>
bool	Container<Element>::remove() {	// deletes current, but does not delete the element
	if(current==NULL) {
		cout<<"CAN'T REMOVE FROM AN EMPTY CONTAINER\n";cout.flush();
		return false;//failed to remove: the list is empty
	}
//-cout<<"Container:"<<this<<": remove:nptrs="<<nptrs<<" ***\n";cout.flush();///DDD
	if(--nptrs==0) {
		GridContainer::pool->put(current);
		current=first=last=NULL;
		return true;
	}
	if(current==last)last=last->prev;
	if(current==first)first=first->next;
	Ptr<Element> *old=current;
	current=current->prev;
	GridContainer::pool->put(old);
//-if(first->next==NULL){cout<<"Container::remove:GOTCHA!!!!!\n";cout.flush();exit(1);}///DDD
//-if(current==GridContainer::pool->hook->next||current->next==GridContainer::pool->hook->next){cout<<"Container::remove:hook->next=current="<<current<<endl;cout.flush();exit(1);}///DDD
	return true;
}
template <class Element>
bool	Container<Element>::remove(Ptr<Element> *ptr) {	// deletes current
//-cout<<"Container:"<<this<<": remove:nptrs="<<nptrs<<", ptr="<<ptr<<endl;cout.flush();///DDD
	if(current==NULL||ptr==NULL) {
		cout<<"Failed to remove from the list: the list is empty\n";cout.flush();
		return false;
	}
//-if(GridContainer::pool->hook->next==ptr){cout<<"Container::remove("<<ptr<<"):hook->next=ptr="<<ptr<<endl;cout.flush();exit(1);}///DDD
	if(nptrs==1){
		if(ptr==first) {
			GridContainer::pool->put(ptr);
			first=last=current=NULL;
			nptrs=0;
			return true;
		} else {
			cout<<"CAN'T REMOVE ptr="<<ptr<<" FROM THE CONTAINER\n";cout.flush();
			return false;//or dump and abort
		}
	}
	if(ptr==current) current=current->next;
	if(ptr==last)last=last->prev;
	if(ptr==first)first=first->next;
	GridContainer::pool->put(ptr);
	nptrs--;
//+	current=ptr; //fast but unsafe
//+	remove();
	return true;
}
//-template <class Element>
//-int	Container<Element>::erase() {	// deletes current
//-	if(current==NULL) return 0;//failed to erase: the list is empty
//-	if(current->next==current) {	
//-		//-delete current->element;
//-		first=last=NULL;
//-		delete current;
//-		current=NULL;
//-		nptrs=0;
//-		return 1;
//-	}
//-	current->next->prev=current->prev;
//-	current->prev->next=current->next;
//-	if(last==current)last=last->prev;
//-	Ptr<Element> *old=current;
//-	current=current->prev;
//-//-	delete old->element;
//-	delete old;
//-	nptrs--;
//-	first=last->next;
//-	return 1;
//-}
//-template <class Element>
//-bool	Container<Element>::erase(Element *element) {
//-	if(!locate(element))return false;
//-	erase();
//-	return true;
//-}
//-template <class Element>
//-bool Container<Element>::locate(Element *element) {
//-	if(current==NULL) return false;
//-	current=first;
//-	do {
//-		if(&current->element==element) return true;
//-		current=current->next;
//-	}	while(current!=first);
//-	return false;
//-}
//-template <class Element>
//-bool Container<Element>::locate(Ptr<Element> *ptr) {
//-	if(current==NULL) return false;
//-	current=first;
//-	do {
//-		if(current==ptr) return true;
//-		current=current->next;
//-	}	while(current!=first);
//-	return false;
//-}

