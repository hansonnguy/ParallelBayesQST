function QqubitBures_v5
% This function performs Bayesian state estimation for multi-qubit
% experiments simulated by simBuresQubits2.m.

%v5 features: parallelism and returns every sample

% JML/HHN
% 2024.06.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
clear all; close all;
%% LOOP SETTINGS
th = 12;                                   % Thinning value.
ch = 0;                                 % Chain value

fileNum = '002';   % Saved file name
A = load('2024-08-09ibmqQubit_Q=4L=1P=400','counts','rhoVecG','Q','P');
counts = A.counts;
rhoVecG = A.rhoVecG;                    % Ground truth.
Q = A.Q;
P = A.P;
clear A;          % Clear up memory space (might not be necessary).

%% INPUTS
numSamp = 2^10;      % Number of samples to obtain from MH.
D = 2^Q;            % Hilbert space dimension.
numCores = 6;      % Core value                                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAP FROM DENSITY MATRIX TO PROBABILITIES:
% Eigenkets of Pauli basis (X,Y,Z):
ket(:,:,1) = [ 1/sqrt(2)*[1 1]; 1/sqrt(2)*[1 -1] ].';       % X
ket(:,:,2) = [ 1/sqrt(2)*[1  1i]; 1/sqrt(2)*[1 -1i] ].';      % Y
ket(:,:,3) = [ [1 0]; [0 1] ].';                            % Z

% Obtain all basis combinations & state labels:
indM = zeros(3^Q,Q);         % Matrix of indices for all possible Pauli combinations.
indS = zeros(2^Q,Q);         % Matrix of indices for all possible outcomes for a given measurement.
for q=1:Q
    indM(:,q) = reshape(repmat((1:1:3).',[3^(q-1) 3^(Q-q)]).',[3^Q 1]);
    indS(:,q) = reshape(repmat((1:1:2).',[2^(q-1) 2^(Q-q)]).',[2^Q 1]);
end

% Form measurement matrix:
M = [];
for m=1:3^Q                 % Cycle through measurements.
    for s=1:2^Q             % Cycle through outcomes.
        KET = ket(:,indS(s,1),indM(m,1));       % Ket for first qubit's state.
        for q=2:Q
            KET = kron(KET,ket(:,indS(s,q),indM(m,q)));     % Kronecker product of all eigenstates for given (m,s)
        end
        proj = kron(conj(KET),KET);             % rho projector.
        M = [M; proj.'];
    end
end


%% INFERENCE LOOP
[~,numStates] = size(rhoVecG);          % Find # of states.
% Fb = zeros(numStates,1);                % Initialize fidelity.
% froB = zeros(numStates,1);              % Initialize Frobenius error.
% samplerTime = zeros(numStates,1);
% rhoB = zeros(D,D,numStates);
ftnHandle = @QqubitNest;

[rho, FbofBayeToGT, thinTime] = parallel_chain(ch, th, numStates);

% WRITING TO FILE
Today = date;
FileName = ['ibmq_QqubitBures_v5_' datestr(Today,'yyyy') datestr(Today,'mm') ...
    datestr(Today,'dd') '_Q=' num2str(Q) '_ch=' num2str(ch) '_th=' num2str(th) '_numSamp=' num2str(numSamp) '_' fileNum];
save(FileName,'rho', 'FbofBayeToGT','thinTime', 'th','ch', 'numSamp', 'numStates','Q')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUB-ROUTINES
    function  [Fb,froB,rhoB, rhoSamps, samplerTime] = QqubitNest(counts,rhoG,THIN)
        %   Fb = fidelity of Bayesian mean wrt ground truth.
        %   froB = squared Frobenius error of Bayesian mean wrt ground truth.
        %   rhoB = Bayesian mean density matrix.
        %   rhoSamps = each sample
        %   counts = ROW vector of counts for all projections.
        %   rhoG = ground truth density matrix.
        %   THIN = thinning value for MCMC.

        
        %% INITIAL VALUES FOR SAMPLER
        rng('shuffle')
        param0 = randn([1 4*D^2]);
        beta = 0.1;             % Initial parameter for stepsize.
        
        %% MH SAMPLER
        samp = zeros(numSamp,4*D^2);    % Initializing samples.
        acc = 0;                            % Counter of acceptances.

        r = 1.1;        % Acceptance rate update factor.
        Mb = 500;       % Update stepsize after this many points.
        
        % Initial point:
        x = param0;
        rhoX = paramToRhoCol(x);              % We draw a random vector, but make sure it is in the space of density matrices
        logX = counts*log(M*rhoX);            % Likelihood.
        
        tic;
        for j=1:numSamp*THIN
            % Proposed updated parameters (all are normal):
            y = sqrt(1-beta^2)*x + beta*randn([1 4*D^2]);
            
            rhoY = paramToRhoCol(y);
            logY = counts*log(M*rhoY);
            
            if log(rand) < logY - logX
                x = y;      % Accept new point.
                logX = logY;
                acc = acc+1;
            end
            
            if mod(j,Mb)==0         % Stepsize adaptation
                rat = acc/Mb;        % Estimate acceptance probability, and keep near 0.234
                if rat > 0.6
                    beta = beta*r;
                    if beta > 1     % Cap beta at unity.
                        beta = 1;
                    end
                elseif rat < 0.2
                    beta = beta/r;
                end
                acc=0;
            end
            
            if mod(j,THIN) == 0
                samp(j/THIN,:) = x;     % Store samples.
            end
        end
        samplerTime = toc;
        
        %% OUTPUT
        rhoB = zeros(D,D);
        rhoSamps = zeros(D,D,numSamp);
        for n=1:numSamp
            rhoAsVec = paramToRhoCol(samp(n,:));
            rhoSamp = reshape(rhoAsVec,[D D]).';    % Watch map convention here.
            rhoSamps(:,:,n) = rhoSamp;
            rhoB = rhoB + 1/numSamp*rhoSamp;        % Bayesian mean.
        end
        Fb = real(trace(sqrtm(sqrtm(rhoG)*rhoB*sqrtm(rhoG)))^2);
        froB = real(trace((rhoB-rhoG)'*(rhoB-rhoG)));
    end

    % Converts parameter set (row vector) into density matrix (expressed as a column vector).
    function z = paramToRhoCol(par)
        G = reshape(par(1:D^2),[D D]) ...
            + 1i*reshape(par(D^2+1:2*D^2),[D D]);        % Ginibre draw 1.
        H = reshape(par(2*D^2+1:3*D^2),[D D]) ...
            + 1i*reshape(par(3*D^2+1:4*D^2),[D D]);      % Ginibre draw 2.
        
        % Find unitary from H (from QETLAB)
        [T,R] = qr(H);
        R = sign(diag(R));
        U = T.*R.';
        W = (eye(D) + U)*G;         % Bures construction.
        
        rho = W*W'/trace(W*W');             % Density matrix.
        z = reshape(rho.',[],1);            % Convert to column, in order of rho(1,1), rho(1,2), rho(1,3),...,etc.
    end



    function [rhoChainSamps, FbofBayeToGT, thinTime] = parallel_chain(ch, th, numStates)

        totChains = 2^ch;

        rhoChainSamps = zeros(D, D, numSamp, totChains, th + 1, ...
            numStates);
        FbofBayeToGT = zeros(1, numStates);  
        thinTime = zeros(1, th + 1); 

        delete(gcp('nocreate'));

        parpool(min(totChains,numCores));

        % Parallel MCMC within each state, thinning, and chain
        for state = 1:numStates 
                
            rhoG = reshape(rhoVecG(:,state),[D D]).';       % Ground truth as matrix.
            rhoChain = zeros(D,D, totChains, th + 1); % store the Bayesian mean of each chain for specified thinning value
            for thin = 1:(th + 1)
              
                tic;
                parfor chain = 1:(totChains)
                    
                    % Markov chains in parallel for each thin
                    [~,~,rhoB,rhoSamps,~] ... 
                    = feval(ftnHandle,counts(state,:),rhoG,2^(thin-1));
                    
                    rhoChain(:,:, chain, thin, state) = rhoB; 
                    rhoChainSamps(:, :, :, chain, thin, state) = rhoSamps; 
                end

                thinTime(1, thin) = thinTime(1, thin) + toc/numStates;
                fprintf(['THIN: ' num2str(thin  - 1) ' of ' num2str(th) ' completed \n']);
            end

                rhoBT = sum(rhoChain(:, :, :, th + 1), 3) / (totChains); % Calculate the Bayesian Truth estimate using all chains at the highest thinning
                FbofBayeToGT(1, state) = real(trace(sqrtm(sqrtm(rhoG)*rhoBT*sqrtm(rhoG)))^2)
                (['<strong> STATE: ' num2str(state) ' of ' num2str(numStates) ' completed </strong> \n']);
        end
    end
end