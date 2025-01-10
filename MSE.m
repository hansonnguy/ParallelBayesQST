function MSE
% This function returns MSE with mixed or unmixed chains for multi-qubit
% tomography by QqubitBures_v5_fast_parallel.m

%% LOOP SETTINGS
A = load('QqubitBures_v5_fast_parallel_20240821_Q=4_ch=10_th=12_numSamp=1024_001', 'th', 'ch', 'numSamp', 'numStates', 'rho', 'Q');
B = load('2024-08-09ibmqQubit_Q=4L=1P=400', 'rhoVecG');

th = A.th;
ch = A.ch;

numSamp = A.numSamp;
numStates = A.numStates;
rho = A.rho;
Q = A.Q;

D = 2^Q;
rhoVecG = B.rhoVecG;
rhoG = reshape(rhoVecG(:,1),[D D]).'; % for future reference, use rhoVecG(:,state)


meanSquareError = zeros(ch + 1, th + 1, numStates);
totChains = 2^ch;

%% MIXED OR UNMIXED
mixed = true;
%% rhoBT or rhoG
useRhoG = false;
%% INFERENCE LOOP
for state = 1:numStates
    rhoBT = sum(rho(:, :, :, th + 1, state), 3) / (totChains); % Calculate the Bayesian Truth estimate using all chains at the highest thinning
    for thin = 1:(th + 1) % add 1 to include thin = 0   

        for chainIndex = 1:(ch + 1)
            R = 2^(chainIndex-1);
            % randomize chain selections or keep them linear
            if mixed == true
                r = randsample((2^ch), R);
            else
                r = 1:R;
            end
            meanRhoB = sum(rho(:, :, r, thin, state), 3) / R; 
            if useRhoG == true            
                meanSquareError(chainIndex, thin, state) = (real(trace((meanRhoB-rhoG)'*(meanRhoB-rhoG))));    
            else
                meanSquareError(chainIndex, thin, state) = (real(trace((meanRhoB-rhoBT)'*(meanRhoB-rhoBT))));
            end
        end
    end
end

%Make sure last MSE is 0 (BT - BT)
meanSquareError(ch + 1, th + 1, numStates) = 0;

fileTag1 = 'mixed_chains';   
fileTag2 = 'GT';   
%WRITING TO FILE
Today = date;
FileName = ['ibmq_QqubitBuresv5_MSE_' datestr(Today,'yyyy') datestr(Today,'mm') ...
    datestr(Today,'dd') '_Q=' num2str(Q) '_ch=' num2str(ch) '_th=' num2str(th) '_numSamp=' num2str(numSamp) '_' fileTag1 '=' num2str(mixed) '_' fileTag2 '=' num2str(useRhoG)];
save(FileName,'meanSquareError', 'th','ch', 'numSamp', 'numStates', 'Q');
