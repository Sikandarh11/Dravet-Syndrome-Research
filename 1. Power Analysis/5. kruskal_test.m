function kruskal_test
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Get all EDF files
    controlFiles = dir(fullfile(controlFolderPath, '*.edf'));
    dsFiles = dir(fullfile(dsFolderPath, '*.edf'));

    % Relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGCz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Process files and extract power values for the relevant electrodes
    fprintf('Processing Control EEG files...\n');
    controlData = process_files(controlFiles, controlFolderPath, relevantElectrodes);

    fprintf('Processing DS EEG files...\n');
    dsData = process_files(dsFiles, dsFolderPath, relevantElectrodes);

    % Perform statistical tests (e.g., T-Test, Kruskal-Wallis)
    perform_kruskal_test(controlData, dsData);
end

%% Helper Function: Process files and extract relevant electrode data
function data = process_files(fileList, folderPath, relevantElectrodes)
    data = [];  % Initialize storage
    for k = 1:length(fileList)
        try
            % Read EDF file
            edfFile = fullfile(folderPath, fileList(k).name);
            eegData = edfread(edfFile);

            % Extract data for relevant electrodes
            electrodeData = [];
            for electrode = relevantElectrodes
                if ismember(electrode, eegData.Properties.VariableNames)
                    signal = eegData{:, electrode};
                    signal = cell2mat(signal);  % Convert to numeric array
                    electrodeData = [electrodeData, signal(~isnan(signal))]; %#ok<AGROW>
                end
            end

            % Store average power for all relevant electrodes
            avgPower = mean(electrodeData, 2);
            data = [data; avgPower]; %#ok<AGROW>

            fprintf('Processed %s with %d relevant signals.\n', ...
                    fileList(k).name, size(electrodeData, 2));
        catch ME
            fprintf('Error processing %s: %s\n', fileList(k).name, ME.message);
        end
    end
end

%% Helper Function: Perform Kruskal-Wallis Test
function perform_kruskal_test(controlData, dsData)
    fprintf('Performing Kruskal-Wallis Test...\n');
    try
        [p, tbl, stats] = kruskalwallis([controlData; dsData], ...
                                         [ones(size(controlData)); 2*ones(size(dsData))], 'off');
        fprintf('p-value: %.4f\n', p);
        disp(tbl);
        disp(stats);
    catch ME
        fprintf('Error in Kruskal-Wallis Test: %s\n', ME.message);
    end
end
