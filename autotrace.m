function autotrace
% This function returns an autotrace matrix for multi-qubit
% tomography by QqubitBures_v5.m

clear all; close all;

A = load('QqubitBures_v5_20240827_Q=1_ch=0_th=12_numSamp=1024_001', 'rho', 'th','ch', 'numSamp', 'Q');

rho     = A.rho;
th      = A.th;
ch      = A.ch;
numSamp = A.numSamp;
Q       = A.Q;

lmax = 200;
acf = zeros(th+1, lmax + 1);
Neff = zeros(th+1, 1);
N = zeros(th+1, 1);
IACT = zeros(th+1, 1);
MSE = zeros(th+1, 1);
rhoMean = sum(rho(:, :, :, th + 1), 3)/numSamp;

for TH = 1:(th+1)
    for l = 0:lmax
        for L = 1:(numSamp-lmax)
        acf(TH, l+1) = acf(TH, l+1) + (real(trace((rho(:, :, L, 1, TH, 1)-rhoMean)'*(rho(:, :, L+l, 1, TH, 1)-rhoMean)))); %ACF function
        end
    end
    acf(TH, :) = acf(TH, :)/max(acf(TH, :));

    % Calculate Neff
    IACT(TH) = 1 + 2*sum(acf(TH, 1:end));  % sum from lag 0 to lmax
    Neff(TH) = numSamp*2^(TH)/IACT(TH);  % Neff calculation
    N(TH) = numSamp*2^(TH);

    meanRhoB = sum(rho(:, :, :, 1, TH, 1), 3) / (numSamp);
    MSE(TH) = (real(trace((meanRhoB-rhoMean)'*(meanRhoB-rhoMean))));
end
disp('IACT for thin = 12: ')
disp(IACT(end))
%% MSE Calculations for a single chain

% Plot MSE vs Neff
% figure;
% plot(Neff, MSE, '-o');
% xlabel('Effective Sample Size (Neff)');
% ylabel('Mean Squared Error (MSE)');
% yscale log
% xscale log
% title('MSE vs Neff');
% grid on;

Today = date;
FileName = ['sim_QqubitBuresv5_ACF_' datestr(Today,'yyyy') datestr(Today,'mm') ...
    datestr(Today,'dd') 'Q=' num2str(Q) '_ch=' num2str(ch) '_th=' num2str(th) '_numSamp=' num2str(numSamp)];
% save(FileName, 'lmax', 'acf', 'th', 'numSamp', 'Q', 'IACT');
