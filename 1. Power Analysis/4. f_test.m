function f_test
    % Helper function for F-Test
    function [h, pValue, F, var1, var2] = ftest2(group1, group2)
        var1 = var(group1);
        var2 = var(group2);
        F = var1 / var2;
        n1 = length(group1) - 1;
        n2 = length(group2) - 1;
        pValue = 2 * (1 - fcdf(max(F, 1/F), n1, n2));
        h = pValue < 0.05;
    end

    % Define directories
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Get relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Process files for control and DS groups
    fprintf('Processing Control EEG files...\n');
    controlGroupData = process_files(controlFolderPath, relevantElectrodes);

    fprintf('Processing DS EEG files...\n');
    dsGroupData = process_files(dsFolderPath, relevantElectrodes);

    % Perform F-Test
    fprintf('Performing F-Test...\n');
    [h, pValue, F, varControl, varDS] = ftest2(controlGroupData, dsGroupData);

    % Display Results
    fprintf('F-Test Results:\n');
    fprintf('Hypothesis result (h): %d\n', h);
    fprintf('p-value: %.4f\n', pValue);
    fprintf('F-statistic: %.4f\n', F);
    fprintf('Control variance: %.4f\n', varControl);
    fprintf('DS variance: %.4f\n', varDS);
end

%% Helper Function to Process Files
function groupData = process_files(folderPath, relevantElectrodes)
    files = dir(fullfile(folderPath, '*.edf'));
    groupData = [];
    for k = 1:length(files)
        try
            edfFile = fullfile(folderPath, files(k).name);
            data = edfread(edfFile);
            electrodeData = [];
            for electrode = relevantElectrodes
                if ismember(electrode, data.Properties.VariableNames)
                    signal = cell2mat(data{:, electrode});
                    electrodeData = [electrodeData, signal(~isnan(signal))];
                end
            end
            avgPower = mean(electrodeData, 2);
            groupData = [groupData; avgPower];
            fprintf('Processed %s with %d electrodes.\n', files(k).name, size(electrodeData, 2));
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end
end
