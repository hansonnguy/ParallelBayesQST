# ParallelBayesQST
This repository contains MATLAB source code that performs parallelized Bayesian quantum state tomography through an unorthodox pooling method of pCN Metropolis-Hastings sampling chains. The code is utilized in [], where the theory is more rigorously described. The codes provided are designed to replicate figures in the reference. Input parameters will allow the user to simulate quantum experiments and perform parallelized quantum state tomography. 

The codes are separated by the simulated data set and the experimental data set. The simulated data can be entirely reproduced by inputting the settings produced below. The experimental data is obtained by QISKIT and reformatted to fit the conventions presented in the MATLAB code. The input parameters and simulation settings used in the paper are further elaborated here:   

Fig. 2: Run froB.m with lmax = 200 with:
- dataFileName = 'ParallelQqubitBures__Q=1_chain=1_th=12_numSamp=1024_001';
- dataFileName = 'ParallelQqubitBures__Q=2_chain=1_th=12_numSamp=1024_001';
- dataFileName = 'ParallelQqubitBures__Q=3_chain=1_th=12_numSamp=1024_001';
- dataFileName = 'ParallelQqubitBures__Q=4_chain=1_th=12_numSamp=1024_001';

Fig. 3: Run squaredFroB.m:
- dataFileName = 'ParallelQqubitBures__Q=1_chain=1024_th=12_numSamp=1024_001’; dataFileName2 = 'acfData_Q=1_th=12_numSamp=1024';
- dataFileName = 'ParallelQqubitBures__Q=2_chain=1024_th=12_numSamp=1024_001'; dataFileName2 = 'acfData_Q=2_th=12_numSamp=1024';
- dataFileName = 'ParallelQqubitBures__Q=3_chain=1024_th=12_numSamp=1024_001'; dataFileName2 = 'acfData_Q=3_th=12_numSamp=1024';
- dataFileName = 'ParallelQqubitBures__Q=4_chain=1024_th=12_numSamp=1024_001'; dataFileName2 = 'acfData_Q=4_th=12_numSamp=1024';

Fig. 4: Run fbTime with qMax = 4; fbPrecision = 6;

Fig. 5: in the Experimental Data set, run bar3(rhoB) for:
‘ParallelQqubitBures_Q=1_chain=1024_th=12_numSamp=1024_001’
‘ParallelQqubitBures_Q=2_chain=1024_th=12_numSamp=1024_001’
‘ParallelQqubitBures_Q=3_chain=1024_th=12_numSamp=1024_001’
‘ParallelQqubitBures_Q=4_chain=1024_th=12_numSamp=1024_001’

Please contact Hanson Nguyen at nhanson7@asu.edu with any comments or questions. Thanks for stopping by!
