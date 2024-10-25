function t_test
    % Helper function for Independent-Samples T-Test
    function [h, pValue, ci, stats] = independent_ttest(group1, group2)
        [h, pValue, ci, stats] = ttest2(group1, group2, 'Vartype', 'unequal');
    end

    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Get the list of EDF files in both directories
    controlFiles = dir(fullfile(controlFolderPath, '*.edf'));
    dsFiles = dir(fullfile(dsFolderPath, '*.edf'));

    % Define the relevant electrodes from the paper
    electrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                  "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                  "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Extract and process signal data from both groups
    fprintf('Processing Control EEG files...\n');
    controlData = extract_electrode_data(controlFiles, controlFolderPath, electrodes);

    fprintf('Processing DS EEG files...\n');
    dsData = extract_electrode_data(dsFiles, dsFolderPath, electrodes);

    % Perform the Independent-Samples T-Test
    fprintf('Performing Independent-Samples T-Test...\n');
    try
        [h, pValue, ci, stats] = independent_ttest(controlData, dsData);

        % Display the results
        fprintf('T-Test Results (Control vs DS cases):\n');
        fprintf('Hypothesis test result (h): %d\n', h);
        fprintf('p-value: %.4f\n', pValue);
        fprintf('95%% Confidence Interval: [%.4f, %.4f]\n', ci(1), ci(2));
        fprintf('Degrees of Freedom: %.4f\n', stats.df);
    catch ME
        fprintf('Error during T-Test: %s\n', ME.message);
    end
end

%% Helper Function: Extract Electrode Data from EDF Files
function groupData = extract_electrode_data(fileList, folderPath, electrodes)
    groupData = [];
    for k = 1:length(fileList)
        try
            % Load EDF file and extract data
            edfFile = fullfile(folderPath, fileList(k).name);
            data = edfread(edfFile);

            % Extract signals from relevant electrodes
            signalData = [];
            for e = 1:length(electrodes)
                if ismember(electrodes(e), data.Properties.VariableNames)
                    electrodeData = data.(electrodes(e));
                    numericData = cell2mat(electrodeData);
                    numericData = numericData(~isnan(numericData));
                    signalData = [signalData; mean(numericData)]; %#ok<AGROW>
                else
                    fprintf('Warning: %s not found in %s\n', electrodes(e), fileList(k).name);
                end
            end

            % Store the mean signal value for this file
            groupData = [groupData; mean(signalData)]; %#ok<AGROW>
            fprintf('Processed %s with %d relevant signals.\n', fileList(k).name, length(signalData));
        catch ME
            fprintf('Error processing %s: %s\n', fileList(k).name, ME.message);
        end
    end
end
