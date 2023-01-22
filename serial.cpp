#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <vector>
#include <valarray>
#include <iomanip>
#include <time.h>

using namespace std;

const double PI = 3.141592653589793238460;
 
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;


// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(double *re, double *im, int N)//double re[], double im[]
{
    if (N <= 1) return;
 
    // divide
    int size = N/2*sizeof(double);
    double *even_real;
    even_real = (double*)malloc(size);
    double *even_im;
    even_im = (double*)malloc(size);
    double *odd_real;
    odd_real = (double*)malloc(size);
    double *odd_im;
    odd_im = (double*)malloc(size);
    for(int i = 0; i < N; i++){
        if(i%2==0){
            even_real[i/2] = re[i];
            even_im[i/2] = im[i];
        }
        else{
            odd_real[(i-1)/2] = re[i];
            odd_im[(i-1)/2] = im[i];
        }
    }
 
    // conquer
    fft(even_real, even_im, N/2);//even_real, even_im
    fft(odd_real, odd_im, N/2);//odd_real, odd_im
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        double COS = cos(-2 * PI * k / N);
        double SIN = sin(-2 * PI * k / N);
        //Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];//(cos0 + isin0)(x + iy) = (...)
        double t_real = odd_real[k]*COS - odd_im[k]*SIN;
        double t_im = odd_im[k]*COS + odd_real[k]*SIN;
        //double sin0 = math.sin(-2 * PI * k / N)
        //use even_real, even_im, odd_real, odd_im
        //x[k    ] = even[k] + t;
        re[k] = even_real[k] + t_real;
        im[k] = even_im[k] + t_im;
        //x[k+N/2] = even[k] - t;
        re[k+N/2] = even_real[k] - t_real;
        im[k+N/2] = even_im[k] - t_im;
    }
    free(even_real);
    free(even_im);
    free(odd_real);
    free(odd_im);
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
    clock_t start, end;
    start = clock();
    
    int sr=4096;
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


    for(int i = 0; i < 26; i++){

        Complex test[sr];
        for(int j=0;j<sr;j++){
        test[j]=ar[i*sr+j];
        }

        CArray s(test, sr);// have to change.
        
        //convert for fft
        double *re;
        re = (double*)malloc(sr*sizeof(double));
        double *im;
        im = (double*)malloc(sr*sizeof(double));
        for(int i = 0; i < sr; i++){
            re[i] = real(s[i]);
            im[i] = imag(s[i]);
        }
        //forward fft
        fft(re, im, sr);

        //convert back to CArray
        for(int i = 0; i < sr; i++){
            Complex temp(re[i], im[i]);
            s[i] = temp;
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
    outdata.close();
    for(auto a:dft)
    {
     	cout<<a.first<<" "<<a.second<<"\n";
    }
    end = clock();
    double time_taken = double(end-start)/double(CLOCKS_PER_SEC);
    cout << "Time taken: " << time_taken << setprecision(6) << endl;
    return 0;
}
