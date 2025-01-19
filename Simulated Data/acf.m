function acf
% This function returns an autocorrelation function for single-chain, multi-qubit
% tomography by ParallelQqubitBures.m

% HHN
% 2024.06.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
clear all; close all;

%% INPUTS
dataFileName = 'ParallelQqubitBures__Q=2_chain=1_th=12_numSamp=1024_001';
lmax = 200;                 % maximum lag

%% LOOP SETTINGS
dataFile = load(dataFileName, 'rhos', 'rhoB', 'th', 'numStates', 'numSamp', 'Q');

rhoB = dataFile.rhoB;
rhos = dataFile.rhos;
th = dataFile.th;
numSamp = dataFile.numSamp;
numStates = dataFile.numStates;
Q = dataFile.Q;

%% AUTOCORRELATION LOOP
acf = zeros(th+1, lmax + 1, numStates);    % initialize autocorrelated function
IACT = zeros(th+1, 1, numStates);          % initialize integrated autocorrelation time
nEff = zeros(th+1, 1, numStates);          % initialize effective number of samples

for state = 1:numStates
    for TH = 1:(th+1)
        for l = 0:lmax
            for L = 1:(numSamp-lmax)
                acf(TH, l+1, state) = acf(TH, l+1, state) + (real(trace((rhos(:, :, L, 1, TH, state)-rhoB(:, :, state))'*(rhos(:, :, L+l, 1, TH, state)-rhoB(:, :, state))))); % ACF function
            end
        end

    acf(TH, :, state) = acf(TH, :)/max(acf(TH, :));
    IACT(TH, 1, state) = 1 + 2*sum(acf(TH, :, state), 2); 
    nEff(TH, 1, state) = numSamp*2^(TH)/IACT(TH, 1, state); 
    end
end

%% WRITING TO FILE
Today = date;
FileName = ['acfData_' datestr(Today,'yyyy') datestr(Today,'mm') ...
    datestr(Today,'dd') '_Q=' num2str(Q) '_th=' num2str(th) '_numSamp=' num2str(numSamp)];
save(FileName, 'lmax', 'acf', 'th', 'numSamp', 'numStates' ,'Q', 'IACT', 'nEff');


%% PLOT OUTPUT
hold on;    box on;
colorOffSet = 2;
C = colororder(turbo(th+1+2*colorOffSet));
for TH = 1:th+1
    color = C(TH + colorOffSet, :);
    plot(0:lmax, acf(TH,:,1), 'Color', color, DisplayName=sprintf('$2^{%d}$', (TH-1)));
end
xlabel('lag');
ylabel('ACF');
axis([0 lmax -0.1 1]);
lgn = legend('Interpreter', 'latex');
lgn.Title.String = 'Thin Factor $T$';
set(gca,'FontSize',12);