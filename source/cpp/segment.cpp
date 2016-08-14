//CNV-TV Implementation using linked lists 
#include<iostream>
#include<cstdio>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<cfloat>
#include<mpi.h>
#include<cstdlib>
using namespace std;

// structure for node
typedef struct node{
	double x;		// contains the value of yi
	struct node *next;
}node; // structure for each element in the set

// strucuture for set
typedef struct set{
	int c;			// number of elements in the set 
	double xc; 		// average of all yi 
	double x;		// xi in the algorithm 
	double X;		// Xi in the algorithm
	double d;		// d in the algorithm 
	struct node *list;
	struct set *next;
}set;
set *head, *tail;		// start and end pointers of the set 
double lambda=0;			// lambda intialized to 0
double vari, sofs=0.0, avg=0.0;	// variance sum of squares and average
vector<double> sicx;			// contains xic corresponding to the min lambda 

/* functions used in the program */
set* insert_set(double); 	// intial insert into the set 
node* insert_node(double);	// insert into the node
set* join(set*);		// joining set(k) and set(k+1)
void calculate(set*,set*,int);	// for the main calculations
double SIC(set *,int);		// function for calculating SIC
void savesicx(set *);		// function for saving xic values
void printstate(set *);		// prints the current state
double rvalue(double d);	// sgn function
int main(int argc,char **argv)
{
	int rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char intstr[100],outstr[100];
	sprintf(intstr,"%s_inp%d",argv[1],rank);
	sprintf(outstr,"%s_out%d",argv[1],rank);
	ifstream inpfile;
	ofstream outfile;
	int n=1;
	int i,k;
	double j;
//	int nzero=0;
	inpfile.open(intstr);
	outfile.open(outstr);
	inpfile>>j;
	avg+=j;
	sofs+=j*j;
	
	// step 1 of CNV-TV algorithm 	
	head=(set *) malloc(sizeof(set));
	head->c=1;
	head->next=NULL;
	head->list=insert_node(j);	// inserting the list with one element 
	head->xc=head->X=j;
//	cout<<j<<endl;
	tail=head;
	while(inpfile>>j)
	{
		
		
			avg+=j;				// avg 
			sofs+=j*j;				// sum of sqaures 
			tail->next=insert_set(j);		// call insert function 		
			tail=tail->next;	
			n++;			// keep track of tail pointer
		
		
	}
	inpfile.close();
	avg=avg/n;				// calculating average
	sofs=sofs/n;				// diving sum of squares of n
	vari=(sofs-(avg*avg));				// calculating variance
	sicx.resize(n);
//	printstate(head);
///////////////////////////////////////////////////////	
	calculate(head,tail,n);
	i=0;
/*	while(sicx[i]!=-1.11)
	{
		cout<<sicx[i]<<endl;
		i++;
	}*/
	for(i=0;i<n;i++)
		outfile<<sicx[i]<<endl;
	MPI_Finalize();
// end of step 1 of CNV-TV algorithm 	
	return 0;
}

// creating of sets
set* insert_set(double i)
{

	set *temp;
	temp=(set *)malloc(sizeof(set));
	temp->next=NULL;
	temp->c=1;
	temp->list=insert_node(i);		// inserting the list with one element
	temp->xc=temp->X=i;			// intializing variables
	return temp;

}
// creation of a single element in the set 
node* insert_node(double i)
{
	node *temp;
	temp= (node*) malloc(sizeof(node));
	temp->next=NULL;
	temp->x=i;
	return temp;				// return the pointer to the caller function 
}
// function to join two sets k and k+1 , input is address of k 
set* join(set* s)
{
	node *t;
	
	if(s->next)
	{
		t=s->list;
		while(t->next!=NULL)
			t=t->next;		// traverse until end of list in k
		t->next=s->next->list;
		s->c=s->c+s->next->c;			// increase the count
		set *temp;
		node *tem;
		tem=s->list;
		s->xc=0;
		while(tem!=NULL)			// calculatin the value of xi~ ie is avg of all yi
		{
			s->xc+=tem->x;		// xic=signma yi	
			tem=tem->next;
		}
		s->xc=s->xc/s->c;				// dividing with  c
		temp=s->next;
		s->next=s->next->next;
		if(temp)	
			free(temp);				// freeing the space
		if(s->next==NULL)			// changing the tail if necessary 
			tail=s;
	}
	return tail;
}
/* Main calculation function f */
void calculate(set* h,set *t,int n)
{
	int k;					// step number
	int i,j;
	double smin=DBL_MAX,sc;			// smin contains minimum of sj
	set *c,*p;				// temporary pointers
	set *vk;				// vk in the algorithm address of the set k
	double sic;				// sic value 
	double msic=DBL_MAX;			// minimum sic value until now
	double smins[n];			// contains all sminimums and the ptr1 is the count  in the smins array
	set *vks[n];
	int ptr1,ptr2;
	
	for(k=0;h!=t;k++)
	{	
		ptr1=0;
		ptr2=0;
		smin=DBL_MAX;
		/* step 2 calculating di*/
		h->d=rvalue(h->next->X-h->X);
//		cout<<"di/c:"<<h->d<<"/"<<h->c<<" ";
		h->d=h->d/h->c;
		p=h;
		c=h->next;
		while(c->next!=NULL)
		{
			c->d=rvalue(p->X-c->X)+rvalue(c->next->X-c->X);		
//			cout<<c->d<<"/"<<c->c<<" ";
			c->d=c->d/c->c;
			p=c;
			c=c->next;
		}
		t->d=rvalue(p->X-t->X);
//		cout<<t->d<<"/"<<t->c<<endl;
		t->d=t->d/t->c;
		/* end of step 2*/
		
		/* step 3 begins*/
		c=h;
//		cout<<"si:";
		while(c->next!=NULL)
		{
			if (c->next->d-c->d==0)
				sc=DBL_MAX;
			else
				sc= (c->X-c->next->X)/(c->next->d-c->d);		// calculating sj
			smins[ptr1++]=sc;						// storing all the scs for finding all minimum scs
			if(sc<smin)
			{
				smin=sc;
				vk=c;
			}
//			cout<<sc<<" ";
			c=c->next;
		}
//		cout<<endl;

		/*step 3 ends here*/
		
		/* step 4 begins */
		c=h;
		while(c!=NULL)
		{
			c->X=c->X+c->d*smin;		// Xj=Xj+dj*s
			c=c->next;
		}
		/* step 4 ends*/
		
		/* step 5 begins  joining multiple minimums*/
		c=h;
		for(i=0;i<ptr1;i++)
		{
			if(smins[i]==smin)
				vks[ptr2++]=c;
			c=c->next;
		}
		for(i=ptr2-1;i>=0;i--)
		{
			t=join(vks[i]);				// adding and deleting 
		}
//		cout<<"p1 "<<ptr1<<" p2 "<<ptr2<<endl;
		/* step 5 ends*/
	/* SIC step after the joing takes place*/
		c=h;					//JUST
		sic=SIC(c,n);
		if(sic<msic)
		{
			msic=sic;
			c=h;				// JUST
			savesicx(c);	//  saves the xi values corresponding to minimum obtained
		}
//		printf("%d %lf\n",k+1,sic);
	/* SIC calculation and saving ends here*/	
		/* step 6 bigns*/
		c=h;					
		while(c!=NULL)
		{
			c->x=c->X;			// xi= Xj 
			c=c->next;
		}
		lambda= lambda + smin;
		c=h;					//JUST
		/* step 6 ends */
	}
	return ;
}
/* Function that returns the SIC value*/
double SIC(set* h,int n)
{

	double sic;				// sic value
	int m=0;				// m number of sets
	set *temp;				// temporary variable
	double sofd=0;				// sum of (yi-xi~)^2			
	double t;				// temporary variable
	temp=h;
	node *temp1;				// temporary variable for list
	while(temp!=NULL)
	{
		m++;
		temp1=temp->list;
		while(temp1!=NULL)
		{
			t=temp1->x-temp->xc;		// yi-xic
			sofd=sofd+t*t;				// adding to sic all the (yi-xic)^2
			temp1=temp1->next;
		}
		temp=temp->next;
	}
//	cout<<"SEA"<<sofd<<" "<<sofd/vari<<" "<<m<<endl;
	sic=(m)*log(n)+(sofd/vari);			//sic=mlnn+sum of squares of differeces 
//	sic=sofd/vari;	
	return sic;	
}

/* Function for saving the xic values */
void savesicx(set *h)
{
	set *c;			// temporary variable
	int ptr=0;		// ptr 
	c=h;
	node *t;
	while(c!=NULL)
	{
		t=c->list;
		while(t!=NULL)
		{
			sicx[ptr++]=c->xc;			// save xic values
			t=t->next;
		}
		c=c->next;
	}
//	sicx[ptr]=-1.11;
	return;
}

/* function to print the state */
void printstate(set *h)
{
	set *t1;
	node *t;
	t1=h;
//	cout<<"Printing state"<<endl;
	while(t1!=NULL)
	{
	
	//	printf("%lf %d ",t1->xc,t1->c);
		t=t1->list;
		while(t!=NULL)
		{
			printf("%lf\n",t1->x);
			t=t->next;
		
		}
//		printf("\n");
		t1=t1->next;
	}
//	cout<<"State Printed"<<endl;
	return ;
}
double rvalue(double d)
{
	if(d==0.0)
		return d;
	else if(d>0.0)
		return 1.0;
	else 
		return -1.0;
}

