#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <vector>
#include <valarray>
#include <omp.h>

using namespace std;

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

//implement forward real DFT
void rdft(double *re, double *im,Complex *x){

    for(int k=0;k<2048;k++){
        re[k]=0;
        im[k]=0;
    }

    for(int k=0;k<2048;k++){
        #pragma omp parallel for shared(re,x,im)
        for(int i=0;i<2048;i++){
            re[k]+=x[i].real()*cos(2*PI*k*i/2048);
            im[k]+=-x[i].real()*sin(2*PI*k*i/2048);
        }
    }
}


double magnitude(complex<double> p, int N)
{ 
	return 2*sqrt(pow(p.real(),2) + pow(p.imag(),2))/N;
}




/* truncate very small numbers to 0 */
double approx_zero(double d) 
{ 
	if (abs(d) < 0.0000000000001)
		return 0;
	else 
		return d;
}





int main(){
    
    int sr=2048;
    //enter sampling rate 

    double ar[sr*26];
    //26 second song, sr data points per second, or sr hz sampling frequency
    ifstream inFile;
    inFile.open("data.txt");
    int i = 0;
    double value;
    while(inFile >> value){
        ar[i] = value;
        i++;
    }

    //frequency amplitude pair
    pair<int,double> P;
    //vector containing max amplitude with frequency
    vector<pair<int,double>> dft;

    ofstream outdata;
	outdata.open("data2.txt");


   //dft part starts
    for (int i = 0; i < 26; i++)
    {
        Complex test[2048];
        for (int j = 0; j < sr; j++)
        {
            test[j] = ar[i * sr + j];
        }
        CArray s(test,sr);//have to change
        
        //forward fft
        double *re;
        re = (double*)malloc(sr*sizeof(double));
        double *im;
        im = (double*)malloc(sr*sizeof(double));
        for(int i = 0; i < sr; i++){
            re[i] = 0;
            im[i] = 0;
        }
        //forward dft on real values
        rdft(re,im,test);

        //convert back to CArray
        for(int i=0;i<sr;i++){
            Complex temp(re[i],im[i]);
            s[i]=temp;
        }

	  int idx=0;
	    outdata.precision(4);
        outdata<<"iteration: "<< i+1 <<endl;
	    for (idx=0; idx < sr; idx++)
		  outdata << idx << " ";
	    outdata << endl;
	    for (idx=0; idx < sr; idx++)
	    	outdata << (magnitude(s[idx],sr)) << " ";
	    outdata << endl <<endl;
        
        int xk = 0;
        int abs=0.0000000000001;
        for(int j = 1; j <sr/2; j++){
            if(magnitude(s[xk],sr)-(magnitude(s[j],sr))<abs)
                xk=j;
                //Add real part
                //xk += ar[i*1000 + j] * cos(2*M_PI*k*j/1054);
            }
        dft.push_back(make_pair(xk,magnitude(s[xk],sr)));
        //Find 1 highest frequency in each second..   
    }

    return 0;
}
