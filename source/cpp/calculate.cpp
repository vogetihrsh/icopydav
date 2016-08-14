#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

double calculateAverage(vector<float>& v){
	int SIZE,i,j,k;		// normal variables
	int NUM_ITER=6;				// number of iterations 
	double ut,lt;				// upper and lower thresholds
	long double avg=0,prev=0;		// return value avg and previous value of avg
	int csize=0;			// current number of bins under consideration 
	
	
	SIZE = v.size();			// size of the input numbers 
	
	for(j=0;j<SIZE;j++)
		avg+=v[j];
	
	avg = (avg*1.0)/(SIZE*1.0);
	ut = avg * 1.5 ;
	lt = avg * 0.5 ;

	for(i=0;prev!=avg && i<NUM_ITER;i++){
		
		prev = avg;
		avg=0;csize=0;

		for(j=0;j<SIZE;j++){
			if(v[j]>ut || v[j]<lt) ;
			
			else {
				avg+=v[j];
				csize++;
			}
		}
		avg=(avg*1.0)/csize;
		ut = avg*1.5; lt = avg*0.5 ;
	}
	return avg;
}
int main(int argc,char** argv)
{
	int SIZE;
	double avg;
	
	SIZE = atoi(argv[1]);
	
	vector<float> v(SIZE);
	int i;
	for(i=0;i<SIZE;i++)
		cin>>v[i];
	avg = calculateAverage(v);
	cout<<avg<<endl;
	return 0;
}
