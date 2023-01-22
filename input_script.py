import librosa, sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
start_time = time.time()
os.chdir(r'C:\Users\Ricky\Desktop\final_fft')
y, sr = librosa.load('hibou_test.wav', sr=4096)
#plt.plot(np.linspace(1,26,26*1<<14), y)
#plt.show()
np.savetxt('data.txt', y)


#26 second song, 1000 data points per second, or 1000hz sampling frequency

def DFT(x):
    N = len(x)
    n = np.arange(N)
    k = np.linspace(1, N, N).reshape((N,1))
    e = np.exp(-2j * np.pi * k * n / N)

    X = np.dot(e, x)
    return X
    

for i in range(0,26):
    X = DFT(y[i*4096:(i+1)*4096])
    freq = 4096
    '''
    plt.figure(figsize = (8,6))
    plt.plot(np.arange(freq), abs(X))
    plt.title(str(i+1)+'th second')
    plt.xlabel('Freq (Hz)')
    plt.ylabel('DFT Amplitude')
    plt.show()
    '''
print(time.time()-start_time)
