function spearman_correlation
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % List of relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Age data for control and DS patients
    controlAges = containers.Map({'71075', '75872', '87925', '110903', ...
                                   '113566', '121139', '144955', '171679'}, ...
                                  [6.5, 9.1, 7.5, 9.5, 10, 9.2, 11.7, 7.8]);
    dsAges = containers.Map({'23281', '73716', '79426', '111822', ...
                             '138825', '148074', '167741'}, ...
                             [15.11, 9.1, 9.1, 16.5, 6.1, 16.2, 11.3]);

    % Initialize storage for group data and ages
    controlDeltaPower = []; controlThetaPower = []; agesControl = [];
    dsDeltaPower = []; dsThetaPower = []; agesDS = [];

    % Process files for the Control group
    fprintf('Processing Control EEG files...\n');
    [controlDeltaPower, controlThetaPower, agesControl] = process_files(controlFolderPath, relevantElectrodes, controlAges, '.edf');

    % Process files for the DS group
    fprintf('Processing DS EEG files...\n');
    [dsDeltaPower, dsThetaPower, agesDS] = process_files(dsFolderPath, relevantElectrodes, dsAges, '.edf');

    % Perform Spearman Correlation Analysis
    fprintf('Performing Spearman Correlation Analysis...\n');
    try
        % Control group delta and theta
        [rhoDeltaControl, pDeltaControl] = spearman_correlation(controlDeltaPower, agesControl);
        [rhoThetaControl, pThetaControl] = spearman_correlation(controlThetaPower, agesControl);
        fprintf('Control Group Delta Power: rho = %.4f, p = %.4f\n', rhoDeltaControl, pDeltaControl);
        fprintf('Control Group Theta Power: rho = %.4f, p = %.4f\n', rhoThetaControl, pThetaControl);

        % DS group delta and theta
        [rhoDeltaDS, pDeltaDS] = spearman_correlation(dsDeltaPower, agesDS);
        [rhoThetaDS, pThetaDS] = spearman_correlation(dsThetaPower, agesDS);
        fprintf('DS Group Delta Power: rho = %.4f, p = %.4f\n', rhoDeltaDS, pDeltaDS);
        fprintf('DS Group Theta Power: rho = %.4f, p = %.4f\n', rhoThetaDS, pThetaDS);
    catch ME
        fprintf('Error during Spearman Correlation Analysis: %s\n', ME.message);
    end
end

%% Helper Function: Process files and extract δ and θ power
function [deltaPower, thetaPower, ages] = process_files(folderPath, electrodes, ageMap, ext)
    files = dir(fullfile(folderPath, ['*' ext]));
    deltaPower = [];
    thetaPower = [];
    ages = [];

    for k = 1:length(files)
        try
            % Extract patient ID without extension
            [~, fileName, ~] = fileparts(files(k).name);
            if isKey(ageMap, fileName)  % Only process if age data exists
                edfFile = fullfile(folderPath, files(k).name);
                info = edfinfo(edfFile);
                data = edfread(edfFile);

                % Extract data for the relevant electrodes
                electrodeData = [];
                for electrode = electrodes
                    if ismember(electrode, data.Properties.VariableNames)
                        signal = data{:, electrode};
                        signal = cell2mat(signal);
                        electrodeData = [electrodeData, signal(~isnan(signal))]; %#ok<AGROW>
                    else
                        fprintf('Warning: Electrode %s not found in file %s\n', electrode, files(k).name);
                    end
                end

                % Compute δ (delta) and θ (theta) power using band-pass filters
                if ~isempty(electrodeData)
                    % Delta power (0.5 - 4 Hz)
                    deltaPowerData = bandpass(electrodeData, [0.5 4], info.NumSamples(1));
                    deltaPower = [deltaPower; mean(mean(deltaPowerData.^2))]; %#ok<AGROW>

                    % Theta power (4 - 8 Hz)
                    thetaPowerData = bandpass(electrodeData, [4 8], info.NumSamples(1));
                    thetaPower = [thetaPower; mean(mean(thetaPowerData.^2))]; %#ok<AGROW>

                    % Extract age
                    age = ageMap(fileName);
                    ages = [ages; age]; %#ok<AGROW>

                    fprintf('Processed %s: %d signals, Age: %.1f\n', files(k).name, size(electrodeData, 2), age);
                else
                    fprintf('Warning: No valid electrode data found for %s\n', files(k).name);
                end
            else
                fprintf('Warning: Age not found for file %s\n', files(k).name);
            end
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end

    % Additional debugging
    fprintf('Processed %d files: Delta Power Size: %d, Theta Power Size: %d, Ages Size: %d\n', length(files), length(deltaPower), length(thetaPower), length(ages));
end

%% Helper Function: Perform Spearman Correlation
function [rho, pValue] = spearman_correlation(data1, data2)
    [rho, pValue] = corr(data1, data2, 'Type', 'Spearman');
end
