% main_interference_figures_WOLA_CPW_single_file.m
%
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Jan 7, 2023
%


clc
clear
close all


fprintf('Starting main_interference_figures_all_in_one.m ... \n\n')
% global figuresFolder settingsFileName
% Definitions
settingsFileName = 'settingsData.mat';
resultsFolder = 'interference_results';
figuresFolder = [resultsFolder '/figures'];
if ~isdir(resultsFolder)  %#ok
    message = strcat(resultsFolder, ' does not exist.');
    error(message)
end
if ~isdir(figuresFolder)  %#ok
    mkdir(figuresFolder)
end
resultFiles = dir(fullfile(resultsFolder));
% Generate figure:
for fileIndex = 1:length(resultFiles)
    if resultFiles(fileIndex).isdir
        continue
    end
    fprintf('Working on file %s.\n', resultFiles(fileIndex).name)
    resultInfo = split(resultFiles(fileIndex).name, '_');
    typeOFDM = resultInfo{end}(1:end-4);
    fprintf('Contains results for %s-OFDM.\n', typeOFDM)
    dataLoader = load([resultsFolder '/' resultFiles(fileIndex).name]);
    settingsLoader = load(settingsFileName);
    cpVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
    rcVector = dataLoader.rcWindowInterference;
    lineWidth = 1;
    markerSize = 15;
    fontSize = 11;
    switch typeOFDM
        case {'wtx', 'wrx', 'CPwtx', 'CPwrx'}
            fig = figure;
            optVector = dataLoader.optWindowInterference;
            ax = subplot(1, 1, 1);
            semilogy(cpVector, optVector, '*-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            title(strcat(typeOFDM, '-OFDM'))
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
        case {'WOLA', 'CPW'}
            fig = figure;
            % Case A
            % Step 1
            optVectorTxStep1 = dataLoader.optInterferenceCaseAStep1;
            ax = subplot(3, 2, 1);
            semilogy(cpVector, optVectorTxStep1, '*-', 'LineWidth', ...
                lineWidth, 'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            title(strcat(typeOFDM, '-OFDM -- Case A (Tx first)'))
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
            % Step 2
            optVectorTxStep2 = dataLoader.optInterferenceCaseAStep2;
            ax = subplot(3, 2, 3);
            semilogy(cpVector, optVectorTxStep2, '*-', 'LineWidth', ...
                lineWidth, 'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
            % Step 3
            optVectorTxStep3 = dataLoader.optInterferenceCaseAStep3;
            ax = subplot(3, 2, 5);
            semilogy(cpVector, optVectorTxStep3, '*-', 'LineWidth', ...
                lineWidth, 'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
            % Case B
            % Step 1
            optVectorRxStep1 = dataLoader.optInterferenceCaseBStep1;
            ax = subplot(3, 2, 2);
            semilogy(cpVector, optVectorRxStep1, '*-', 'LineWidth', ...
                lineWidth, 'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            title(strcat(typeOFDM, '-OFDM -- Case B (Rx first)'))
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
            % Step 2
            optVectorRxStep2 = dataLoader.optInterferenceCaseBStep2;
            ax = subplot(3, 2, 4);
            semilogy(cpVector, optVectorRxStep2, '*-', 'LineWidth', ...
                lineWidth, 'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
            % Step 3
            optVectorRxStep3 = dataLoader.optInterferenceCaseBStep3;
            ax = subplot(3, 2, 6);
            semilogy(cpVector, optVectorRxStep3, '*-', 'LineWidth', ...
                lineWidth, 'MarkerSize', markerSize), hold on
            semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
                'MarkerSize', markerSize), hold off, grid on
            set(ax, 'FontSize', fontSize)
            set(ax, 'TickLabelInterpreter', 'latex')
            set(ax, 'linewidth', lineWidth)
            set(ax, 'XColor', 'k')
            set(ax, 'YColor', 'k')
            xlabel('Cyclic Prefix Length, $\mu$', 'interpreter', ...
                'latex', 'FontSize', fontSize)
            ylabel('Interference Power', 'interpreter', 'latex', ...
                'FontSize', fontSize)
            lgd = legend('Optimal', 'Raised Cosine');
            lgd.Interpreter = 'latex';
            lgd.FontSize = fontSize;
        otherwise
            fprintf('Skipping %s-OFDM\n', typeOFDM)
    end
    savefig(fig, strcat('interference_power_', typeOFDM, '.fig'), ...
        'compact')
end


% EoF
