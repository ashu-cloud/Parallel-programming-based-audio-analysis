#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <vector>
#include <valarray>

using namespace std;

//const double PI = 3.141592653589793238460;
 
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

__global__ void add(double *re, double *im, int *x, double *odd_real, double *odd_im, double *even_real, double *even_im)
{
    int k = threadIdx.x + blockIdx.x*blockDim.x;//Supposed to count till N/2..
    int N = *x;
    double COS = cos(-2 * 3.141592653589793238460 * k / N);
    double SIN = sin(-2 * 3.141592653589793238460 * k / N);
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

    // taking 8 blocks and N/8 threads

    int num_threads = N/16; //4096/4 = 1024 = max num of threads
    int num_blocks = 8;
    if(N < 16)
    {
        num_threads = N;
        num_blocks = 1;
    }

    double *dev_re, *dev_im, *dev_odd_real, *dev_odd_im, *dev_even_real, *dev_even_im;
    cudaMalloc( (void**)&dev_re, 2*size );
    cudaMalloc( (void**)&dev_im, 2*size );
    cudaMalloc( (void**)&dev_odd_real, size );
    cudaMalloc( (void**)&dev_odd_im, size );
    cudaMalloc( (void**)&dev_even_real, size );
    cudaMalloc( (void**)&dev_even_im, size );
    int *dev_N;
    cudaMalloc( (void**)&dev_N, sizeof(int) );
    
    cudaMemcpy( dev_re, re, 2*size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_im, im, 2*size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_odd_real, odd_real, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_odd_im, odd_im, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_even_im, even_im, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_even_real, even_real, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_N, &N, sizeof(int), cudaMemcpyHostToDevice );


    add<<< num_blocks , num_threads >>>(dev_re, dev_im, dev_N, dev_odd_real, dev_odd_im, dev_even_real, dev_even_im);

    cudaMemcpy( im, dev_im, 2*size, cudaMemcpyDeviceToHost );
    cudaMemcpy( re, dev_re, 2*size, cudaMemcpyDeviceToHost );

    free(even_real);
    free(even_im);
    free(odd_real);
    free(odd_im);

    cudaFree(dev_im);
    cudaFree(dev_re);
    cudaFree(dev_odd_real);
    cudaFree(dev_even_real);
    cudaFree(dev_odd_im);
    cudaFree(dev_even_im);
    cudaFree(dev_N);
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
    
    const int sr=16384;
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

    //vector containing max amplitude with frequency
    vector<pair<int,double>> dft;

    ofstream outdata;
	outdata.open("data2.txt");


    for(i = 0; i < 26; i++){

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
        for(int j = 0; j < sr; j++){
            re[j] = real(s[j]);
            im[j] = imag(s[j]);
        }
        //forward fft
        fft(re, im, sr);

        //convert back to CArray
        for(int j = 0; j < sr; j++){
            Complex temp(re[j], im[j]);
            s[j] = temp;
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
    return 0;
}
