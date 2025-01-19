%% simBuresQubits2.m
% This code is designed to simulate experiments on a random state with an 
% arbitrary number of qubits. It differs from simBuresQubits.m by (i)
% deterministically cycling through Pauli measurements (rather than 
% choosing each setting randomly) and (ii) saving only the total counts.

% JML
% 2021.03.11
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
close all;  clear all;

%% INPUTS
Q = 2;          % Number of qubits in system.
L = 1;         % Number of states.
P = 800;       % Number of meaurements per Pauli setting.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAP FROM DENSITY MATRIX TO PROBABILITIES:
% Eigenkets of Pauli basis (X,Y,Z):
ket(:,:,1) = [ 1/sqrt(2)*[1 1]; 1/sqrt(2)*[1 -1] ].';       % X
ket(:,:,2) = [ 1/sqrt(2)*[1  i]; 1/sqrt(2)*[1 -i] ].';      % Y
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
for m=1:3^Q                 % Cycle through pauli measurements.
    for s=1:2^Q             % Cycle through outcomes.
        KET = ket(:,indS(s,1),indM(m,1));       % Ket for first qubit's state.
        for q=2:Q
            KET = kron(KET,ket(:,indS(s,q),indM(m,q)));     % Kronecker product of all eigenstates for given (m,s)
        end
        proj = kron(conj(KET),KET);             % rho projector.
        M = [M; proj.'];
    end
end 
size(M)
%% CYCLE THROUGH STATES AND BASES
D = 2^Q;                    % Dimensionality of Hilbert space.
% events = zeros(L,6^Q,E);    % 3D array of events: rows are states, columns are projectors, depth is measurements.
% basis = zeros(L,E);         % Listing of bases chosen for each event.
% outcome = zeros(L,E);       % Outcome of each event (number in {1,2,3,4}).
rhoVecG = zeros(D^2,L);     % Ground truth states (density matrix as column vectors.
counts = zeros(L,6^Q);      % Observations for all states.

for m = 1:L
    x = randn([1 4*D^2]);
    
    G = reshape(x(1:D^2),[D D]) ...
        + i*reshape(x(D^2+1:2*D^2),[D D]);        % Ginibre draw 1.
    H = reshape(x(2*D^2+1:3*D^2),[D D]) ...
        + i*reshape(x(3*D^2+1:4*D^2),[D D]);      % Ginibre draw 2.
    
    % Find unitary from H (from QETLAB)
    [T,R] = qr(H);
    R = sign(diag(R));
    U = T.*R.';
    W = (eye(D) + U)*G;         % Bures construction.
    
    rho = W*W'/trace(W*W');                 % Density matrix.
    rhoVecG(:,m) = reshape(rho.',[],1);     % Convert to column, in order of rho(1,1), rho(1,2), rho(1,3),...,etc.

    prob = real(M*rhoVecG(:,m));                  % Compute probabilities.
    
    for n=1:3^Q                 % Cycle through measurement settings.
        outVec = mnrnd(P,prob(2^Q*(n-1)+1:2^Q*n).');    % Simulate P trials.
        counts(m,2^Q*(n-1)+1:2^Q*n) = outVec;          % Save to counts for state m.
    end    
end

%% WRITING TO FILE
% Today = date;
% FileName = ['simBuresQubits2_' datestr(Today,'yyyy') datestr(Today,'mm') ...
%     datestr(Today,'dd') '_Q=' num2str(Q) '_L=' num2str(L) '_P=' num2str(P) '.mat'];
% save(FileName,'rhoVecG','P','counts','Q', '-v7.3')

