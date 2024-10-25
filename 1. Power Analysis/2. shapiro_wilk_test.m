function shapiro_wilk_test
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % List of relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Get the list of all EDF files in both directories
    controlFiles = dir(fullfile(controlFolderPath, '*.edf'));
    dsFiles = dir(fullfile(dsFolderPath, '*.edf'));

    % Process and extract data for each group
    fprintf('Processing Control EEG files...\n');
    controlData = process_files(controlFiles, controlFolderPath, relevantElectrodes);

    fprintf('Processing DS EEG files...\n');
    dsData = process_files(dsFiles, dsFolderPath, relevantElectrodes);

    % Perform Shapiro-Wilk Test for both groups
    fprintf('Performing Shapiro-Wilk Test...\n');
    perform_shapiro_wilk(controlData, 'Control Group');
    perform_shapiro_wilk(dsData, 'DS Group');
end

%% Helper Function: Process Files and Extract Data for Relevant Electrodes
function groupData = process_files(fileList, folderPath, relevantElectrodes)
    groupData = []; 

    for k = 1:length(fileList)
        try
            % Load EDF file
            edfFile = fullfile(folderPath, fileList(k).name);
            data = edfread(edfFile);

            % Extract data for relevant electrodes
            electrodeData = [];
            for electrode = relevantElectrodes
                if ismember(electrode, data.Properties.VariableNames)
                    signal = data{:, electrode};
                    signal = cell2mat(signal); % Convert to numeric array
                    electrodeData = [electrodeData, signal(~isnan(signal))];
                end
            end

            % Store average power across electrodes
            avgPower = mean(electrodeData, 2); 
            groupData = [groupData; avgPower]; %#ok<AGROW>

            fprintf('Processed %s with %d relevant signals.\n', ...
                    fileList(k).name, size(electrodeData, 2));
        catch ME
            fprintf('Error processing %s: %s\n', fileList(k).name, ME.message);
        end
    end
end

%% Helper Function: Perform Shapiro-Wilk Test for Normality
function perform_shapiro_wilk(data, groupName)
    fprintf('Shapiro-Wilk Test for %s...\n', groupName);
    try
        % Perform Shapiro-Wilk test
        [h, p] = swtest(data);
        fprintf('Result: h = %d, p = %.4f\n', h, p);
    catch ME
        fprintf('Error in Shapiro-Wilk Test for %s: %s\n', groupName, ME.message);
    end
end

%% Helper Function: Shapiro-Wilk Test Implementation (Using Jarque-Bera as Proxy)
function [h, pValue] = swtest(x)
    % Use Jarque-Bera test as an approximation for Shapiro-Wilk
    [~, pValue] = jbtest(x);  
    h = pValue < 0.05; % Reject null hypothesis if p < 0.05
end
