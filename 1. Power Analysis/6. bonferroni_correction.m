function bonferroni_correction
    % Helper function to apply Bonferroni correction
    function [adjustedPValues, significant] = bonferroni_correction(pValues, alpha)
        % Adjusts p-values for multiple comparisons
        m = length(pValues); 
        adjustedPValues = min(pValues * m, 1);
        significant = adjustedPValues < alpha;
    end

    % Folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Relevant electrodes from the paper
    electrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                  "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                  "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Read and process data for both groups
    fprintf('Processing Control EEG files...\n');
    controlData = process_group_files(controlFolderPath, electrodes);

    fprintf('Processing DS EEG files...\n');
    dsData = process_group_files(dsFolderPath, electrodes);

    % Define frequency ranges
    frequencyRanges = [1 4; 4 8; 8 13; 13 30]; 

    pValues = []; 

    % Compare frequency ranges
    fprintf('Comparing frequency ranges...\n');
    for i = 1:size(frequencyRanges, 1)
        fprintf('Testing frequency range: %d-%d Hz\n', frequencyRanges(i, 1), frequencyRanges(i, 2));

        % Extract frequency-specific data
        controlFreqData = extract_frequency_data(controlData, frequencyRanges(i, :));
        dsFreqData = extract_frequency_data(dsData, frequencyRanges(i, :));

        % Perform T-Test
        [~, pValue] = ttest2(controlFreqData, dsFreqData, 'Vartype', 'unequal');
        pValues = [pValues; pValue];

        fprintf('p-value for range %d-%d Hz: %.4f\n', frequencyRanges(i, 1), frequencyRanges(i, 2), pValue);
    end

    % Apply Bonferroni correction
    fprintf('Applying Bonferroni Correction...\n');
    alpha = 0.05;
    [adjustedPValues, significant] = bonferroni_correction(pValues, alpha);

    % Display results
    fprintf('Adjusted p-values:\n');
    disp(adjustedPValues);
    fprintf('Significance Status (1 = significant, 0 = not significant):\n');
    disp(significant);
end

%% Helper Function: Process Group Files for Relevant Electrodes
function data = process_group_files(folderPath, electrodes)
    files = dir(fullfile(folderPath, '*.edf'));
    data = [];
    for k = 1:length(files)
        try
            edfFile = fullfile(folderPath, files(k).name);
            eegData = edfread(edfFile);

            % Extract relevant electrode data
            electrodeData = [];
            for electrode = electrodes
                if ismember(electrode, eegData.Properties.VariableNames)
                    signal = eegData{:, electrode};
                    signal = cell2mat(signal); 
                    electrodeData = [electrodeData, signal(~isnan(signal))];
                end
            end
            data = [data; mean(electrodeData, 2)];
            fprintf('Processed %s with %d relevant signals.\n', files(k).name, size(electrodeData, 2));
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end
end

%% Helper Function: Extract Frequency Data
function freqData = extract_frequency_data(data, freqRange)
    % Simulate data filtering based on frequency range (replace with real logic)
    freqData = data(mod(1:length(data), 30) >= freqRange(1) & ...
                    mod(1:length(data), 30) < freqRange(2));
end
