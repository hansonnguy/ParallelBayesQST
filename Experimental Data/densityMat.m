function densityMat
% This function returns a density matrix plot from the results of
% ParallelQqubitBures.m

% HHN
% 2024.06.19
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
clear all; close all;

%% INPUTS
dataFileName = 'ParallelQqubitBures_Q=3_chain=1024_th=12_numSamp=1024_001';

%% PLOT DENSITY MATRIX
dataFile = load(dataFileName, 'rhoB', 'Fb');
rhoB = dataFile.rhoB;
Fb = dataFile.Fb;
hfig = figure;
plot_density(rhoB, ['$\mathcal{F} = ' num2str(Fb) '$']);

% plot settings
set(findall(hfig,'-property','FontSize'),'FontSize',12) % adjust fontsize to your document
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')


    function plot_density(rho, plot_title)

        b_real = bar3(real(rho));
        numQubits = log2(length(rho)); 

        % label ticks
        labels = cell(2^numQubits, 1);
        for i = 0:(2^numQubits - 1)
            labels{i+1} = ['$|' dec2bin(i, numQubits) '\rangle$'];  % Binary representation for each basis state
        end

        tickIdx = [];
        tickLabels = {};

        for i = 0:(2^numQubits - 1)
            binaryStr = dec2bin(i, numQubits); 
            if binaryStr(2:end) == repmat('0', 1, numQubits - 1) 
                tickIdx = [tickIdx, i+1]; 
                tickLabels = [tickLabels; labels{i+1}]; 
            end
        end

        set(gca, 'XTick', tickIdx, 'XTickLabel', tickLabels);
        set(gca, 'YTick', tickIdx, 'YTickLabel', tickLabels);

        title(plot_title);
        custom_color = [1, 0.4, 0.4]; 

        % Set the color of each bar to red
        for k = 1:length(b_real)
            b_real(k).FaceColor = custom_color;
        end
    end
end
