function squaredFroB
% This function returns squared frobenius error mulit-chain, multi-qubit
% tomography by ParallelQqubitBures.m

% HHN
% 2024.06.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
clear all; close all;

%% INPUTS
dataFileName = 'ParallelQqubitBures__Q=4_chain=1024_th=12_numSamp=1024_001';
dataFileName2 = 'acfData_Q=4_th=12_numSamp=1024';

%% LOOP SETTINGS
dataFile = load(dataFileName, 'th', 'numChains', 'numSamp', 'numStates', 'rhos', 'rhoB', 'Q');
dataFile2 = load(dataFileName2, 'IACT');

th = dataFile.th;
ch = log2(dataFile.numChains);

numSamp = dataFile.numSamp;
numStates = dataFile.numStates;
numChains = dataFile.numChains;
rhos = dataFile.rhos;
rhoB = dataFile.rhoB;
Q = dataFile.Q;

IACT = dataFile2.IACT;
clear dataFile; clear dataFile2;

%% SQUARED FROBENIUS ERROR LOOP
squaredFroB = zeros(ch + 1, th + 1, numStates);     % initialize matrix to store squared Frobenius error
for state = 1:numStates
    for thin = 1:(th + 1)
        for chainIndex = 1:(ch + 1)

            R = 2^(chainIndex-1);
            r = randsample((2^ch), R);
            meanRhoB = sum(rhos(:, :, r, thin, state), 3) / R; 
            squaredFroB(chainIndex, thin, state) = (real(trace((meanRhoB-rhoB)'*(meanRhoB-rhoB))));
        
        end
    end
end

squaredFroB(ch + 1, th + 1, numStates) = 0;

%% WRITING TO FILE
Today = date;
FileName = ['squaredFroBData_' datestr(Today,'yyyy') datestr(Today,'mm') ...
    datestr(Today,'dd') '_Q=' num2str(Q) '_chain=' num2str(numChains) '_th=' num2str(th) '_numSamp=' num2str(numSamp)];
save(FileName, 'squaredFroB', 'th', 'numChains', 'numSamp', 'numStates', 'Q');

%% PLOT OUTPUT

dataFile = kron(numSamp*2.^((0:(th)))./(IACT.'), 2.^(0:ch));
dataFile2 = reshape(dataFile, (ch+1), th+1);

hold on;    box on;
colorOffSet = 2;
C = colororder(turbo(th+1+2*colorOffSet));
for TH = 1:th+1
    color = C(TH + colorOffSet, :);
    plot(dataFile2(:, TH), squaredFroB(:, TH, 1), '-o', 'MarkerFaceColor', color, 'Color', color, DisplayName=sprintf('$2^{%d}$', (TH-1))); 
end
xlabel('$$N_\mathrm{eff}$$', 'Interpreter', 'latex');
ylabel('$$\epsilon_{f}^2$$', 'Interpreter', 'latex');
lgn = legend('Interpreter', 'latex');
lgn.Title.String = 'Thin Factor $T$';
set(gca,'FontSize', 12);
yscale log
xscale log
