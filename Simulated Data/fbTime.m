function fbTime
% This function returns fidelity matrix for multi-qubit
% tomography by ParallelQqubitBures.m

% HHN
% 2024.06.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
clear all; close all;

%% INPUTS
qMax = 4;
fbPrecision = 6; % order of magnitude of accuracy (e.g., 0.99999)

simFiles = cell(qMax, 1);   % load simulation files
FbFiles = cell(qMax, 1);    % load fidelity files
for Q = 1:qMax
    simFiles{Q} = load(sprintf('ParallelQqubitBures__Q=%d_chain=1024_th=12_numSamp=1024_001', Q), ...
        'rhos', 'Fb', 'avgTime', 'th', 'numChains', 'numSamp', 'numStates', 'Q');
    FbFiles{Q} = load(sprintf('fbData_Q=%d_numChains=1024_th=12_numSamp=1024', Q), ...
        'Fb', 'th', 'numChains', 'numSamp', 'numStates', 'Q');
end

%% LOOP SETTINGS
numChains = simFiles{1}.numChains;
th = simFiles{1}.th;
ch = log2(numChains);

%% FIDELITY AND TIME LOOP
% Extract fidelities and normalize times
fb = zeros(qMax, ch+1, th+1); % All fidelities
fbP = zeros(qMax, th+1); % Parallel chain fidelities
fbS = zeros(qMax, th+1); % Single chain fidelities
time = zeros(qMax, th+1); % Wall clock times

cores = 48;
timeNormFactor = cores / numChains;
for Q = 1:qMax
    fb(Q, :, :) = FbFiles{Q}.Fb;
    fbP(Q, :) = FbFiles{Q}.Fb(end, :);
    fbS(Q, :) = FbFiles{Q}.Fb(1, :);

    time(Q, :) = simFiles{Q}.avgTime * timeNormFactor; % Normalize time based on core usage
end

% Preallocate arrays for minimum time
minTime = zeros(qMax, fbPrecision);

% Find best times for fidelity tolerances
for qInd = 1:qMax
    for fbInd = fbPrecision:-1:1
        fbT = 10^(-fbInd);      % check shortest times for these fidelities
        timeMax = max(max(time));
        for CH = 1:(ch+1)
            for TH = 1:(th+1)
                fidelity = fb(qInd, CH, TH);
                avgTime = time(qInd, TH); % average time across chains
                if fidelity <= fbT && avgTime <= timeMax
                    timeMax = avgTime;
                    if fbInd < fbPrecision
                        minTime(qInd, fbInd) = min(time(qInd, TH), minTime(qInd, fbInd+1)); % in case of less time at higher fidelity
                    else
                        minTime(qInd, fbInd) = time(qInd, TH);
                    end
                end
            end
        end
    end
end

% Plot setup
C = colororder(turbo(qMax + 2));
outerLayout = tiledlayout(2, 1, "TileSpacing", 'compact');

% Top layout for comparing Fb_s and Fb_r
topLayout = tiledlayout(outerLayout, 1, 4, "TileSpacing", 'tight');
for qInd = 1:qMax
    nexttile(topLayout, qInd);
    hold on; box on;

    color = C(qInd+1, :);

    % Plot single and parallel chain fidelities
    plot(time(qInd, :), fbP(qInd, :), '-o', 'MarkerSize', 4, ...
        'MarkerFaceColor', color, 'Color', color, 'LineWidth', 2, 'DisplayName', '$1$');
    plot(time(qInd, :), fbS(qInd, :), ':^', 'MarkerSize', 4, ...
        'MarkerFaceColor', color, 'Color', color, 'LineWidth', 2, 'DisplayName', ['$' num2str(numChains) '$']);
    title(['$Q = ' num2str(qInd) '$'], 'Interpreter','latex');
    % Configure axes
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, ...
        'LineWidth', 1.5);
    if qInd == 1
        ylabel('$1-\mathcal{F}$', 'Interpreter', 'latex');
    else
        set(gca, 'YTickLabel', []);
    end

    xlim([min(min(time)), max(max(time))]);
    xlabel('Wall clock time (s)', 'Interpreter', 'latex');
    ylim([min(min([fbP fbS])), max(max([fbP fbS]))]);
    lgn = legend('location', 'best' ,'Interpreter','latex');
    lgn.Title.String = 'chains $R$';
    hold off;
end

% Bottom layout for scaling plot
nexttile(outerLayout, [1, 1]);
hold on; box on;
C = colororder(turbo(fbPrecision+1));
for fbInd = 1:(fbPrecision-1)
    color = C(fbInd+1, :);
    DisplayName = ['$k = ' num2str(fbInd) '$'];
    plot(1:qMax, minTime(:, fbInd), '-', 'MarkerFaceColor', color, ...
        'MarkerSize', 4, 'Color', color, 'LineWidth', 2, 'DisplayName', DisplayName);
end

% Finalize plot
ylabel('Time (s)', 'Interpreter', 'latex');
xlabel('Qubits', 'Interpreter', 'latex');
set(gca,'xtick',0:qMax)
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
lgn = legend('location', 'best' ,'Interpreter','latex');
lgn.Title.String = '$\mathcal{F} > 1-10^{-k}$';
hold off;

end
