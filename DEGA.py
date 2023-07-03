from operator import itemgetter
import numpy as np
import math
N = 8
Kp = 4
dsnr_db = 2
snr = np.power(10, dsnr_db / 10)


def phi_inv(x: float):
    
    if (x>12):
        return 0.9861 * x - 2.3152
    elif (x<=12 and x>3.5):
        return x*(0.009005 * x + 0.7694) - 0.9507
    elif (x<=3.5 and x>1):
        return x*(0.062883*x + 0.3678)- 0.1627
    else:
        return x*(0.2202*x + 0.06448)

mllr = np.zeros(N, dtype=float)

sigma_sq = 1/(2*Kp/N*np.power(10,dsnr_db/10))
mllr[0] = 2/sigma_sq


for level in range(1, int(np.log2(N)) + 1):

    B = np.power(2, level)

    for j in range(int(B / 2)):
        T = mllr[j]
        mllr[j] = phi_inv(T)
        mllr[int(B / 2 + j)] = 2 * T


mask = [[i, 0.0, 1] for i in range(N)]
 
reliability = mllr

for i in range(N):
    mask[i][1] = reliability[i]
       
mask = sorted(mask, key=itemgetter(1), reverse=False) 
    
for i in range(N-Kp):
    mask[i][2] = 0

mask = sorted(mask, key=itemgetter(0))

print("posições: {}".format(np.array([i[2] for i in mask])))
    

pe = np.zeros(N, dtype=float)
for ii in range(N):

    pe[ii] = 0.5 - 0.5 * math.erf( np.sqrt(mllr[ii])/2 )
