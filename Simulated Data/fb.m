function fb
% This function returns fidelity matrix for multi-qubit
% tomography by ParallelQqubitBures.m

% HHN
% 2024.06.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
clear all; close all;

%% INPUTS
dataFileName = 'ParallelQqubitBures__Q=2_chain=1024_th=12_numSamp=1024_001';

%% LOOP SETTINGS
dataFile = load(dataFileName, 'th', 'numChains', 'numSamp', 'numStates', 'rhos', 'rhoB', 'Q');

th = dataFile.th;
ch = log2(dataFile.numChains);

numSamp = dataFile.numSamp;
numChains = dataFile.numChains;
numStates = dataFile.numStates;
rhoB = dataFile.rhoB;
rhos = dataFile.rhos;
Q = dataFile.Q;
D = 2^Q;

Fb = zeros(ch + 1, th + 1, numStates);

%% FIDELITY LOOP
for state = 1:numStates
    for thin = 1:(th + 1)

        for CH = 1:(ch + 1)
            R = 2^(CH-1);
            r = randsample((2^ch), R); % randomize chain selections
            meanRhoB = sum(rhos(:, :, r, thin, state), 3) / R; 
            Fb(CH, thin, state) = (1 - (real(trace(sqrtm(sqrtm(rhoB(:,:,state))*meanRhoB*sqrtm(rhoB(:,:,state))))^2)));
        end
    end
end

%% WRITING TO FILE
Today = date;
FileName = ['simData_Fb_' datestr(Today,'yyyy') datestr(Today,'mm') ...
    datestr(Today,'dd') '_Q=' num2str(Q) '_numChains=' num2str(numChains) '_th=' num2str(th) '_numSamp=' num2str(numSamp)];
save(FileName,'Fb', 'th','numChains', 'numSamp', 'numStates', 'Q');

