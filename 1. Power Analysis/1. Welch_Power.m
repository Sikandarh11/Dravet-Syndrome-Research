function welch_power
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Get the list of all EDF files in both directories
    controlFiles = dir(fullfile(controlFolderPath, '*.edf'));
    dsFiles = dir(fullfile(dsFolderPath, '*.edf'));

  
    fprintf('Found %d control EEG files.\n', length(controlFiles));
    fprintf('Found %d DS EEG files.\n', length(dsFiles));

    % Initialize arrays to store normalized power spectra
    controlPower = [];
    dsPower = [];

    % Process files for both groups
    fprintf('Processing Control EEG files...\n');
    controlPower = process_group_files(controlFiles, controlFolderPath);

    fprintf('Processing DS EEG files...\n');
    dsPower = process_group_files(dsFiles, dsFolderPath);

    % Save power data to .mat file for later use in statistical tests
    save('power_spectra.mat', 'controlPower', 'dsPower');
end

%% Helper Function: Process Group Files and Extract Power Spectra
function powerData = process_group_files(fileList, folderPath)
    % Define relevant electrodes for analysis
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];
    
    powerData = [];  % Initialize storage for power data

    for k = 1:length(fileList)
        try
            filePath = fullfile(folderPath, fileList(k).name);
            data = edfread(filePath);

            % Extract EEG signal data for relevant electrodes only
            electrodeData = [];
            for electrode = relevantElectrodes
                if ismember(electrode, data.Properties.VariableNames)
                    signal = data{:, electrode};
                    signal = cell2mat(signal);
                    electrodeData = [electrodeData, signal(~isnan(signal))]; %#ok<AGROW>
                else
                    fprintf('Warning: Electrode %s not found in file %s\n', electrode, fileList(k).name);
                end
            end

            % Compute power spectrum using Welch's method if data is available
            if ~isempty(electrodeData)
                [frequencies, power] = compute_welch_power(electrodeData);

                % Normalize power within the 1-80 Hz range
                normalizedPower = normalize_power(frequencies, power);
                powerData = [powerData; normalizedPower]; %#ok<AGROW>

                % Debugging: Output power spectrum status
                fprintf('Processed %s. Normalized Power: %s\n', ...
                    fileList(k).name, mat2str(normalizedPower));
            else
                fprintf('No data for relevant electrodes in file %s.\n', fileList(k).name);
            end
        catch ME
            fprintf('Error processing %s: %s\n', fileList(k).name, ME.message);
        end
    end
end

%% Helper Function: Compute Power Spectrum Using Welch’s Method
function [frequencies, power] = compute_welch_power(signal)
    % Parameters for Welch's method
    windowLength = 2 * 250;  % 2-second Hann window (assuming 250 Hz sample rate)
    overlap = 0.5 * windowLength;  % 50% overlap

    % Apply Welch's method to compute power spectrum
    [power, frequencies] = pwelch(signal, hann(windowLength), overlap, [], 250);

    % Debugging: Output frequency range
    fprintf('Computed power spectrum. Frequency range: %.2f-%.2f Hz\n', ...
        min(frequencies), max(frequencies));
end

%% Helper Function: Normalize Power Spectra
function normalizedPower = normalize_power(frequencies, power)
    % Frequency ranges of interest
    bands = struct('delta', [1 4], 'theta', [4 7], 'alpha', [8 12], ...
                   'spindle', [12 16], 'beta', [13 30], 'gamma', [30 80]);

    % Initialize array to store normalized power for each band
    normalizedPower = zeros(1, length(fieldnames(bands)));

    % Total power from 1 to 80 Hz
    validIndices = frequencies >= 1 & frequencies <= 80;
    totalPower = sum(power(validIndices));

    % Compute and normalize power for each frequency band
    bandNames = fieldnames(bands);
    for i = 1:length(bandNames)
        band = bands.(bandNames{i});
        bandIndices = frequencies >= band(1) & frequencies <= band(2);
        bandPower = sum(power(bandIndices));
        normalizedPower(i) = log10(bandPower / totalPower);  % Log base 10 transform
    end

    % Debugging: Output normalized power for each band
    fprintf('Normalized power (log-transformed): %s\n', mat2str(normalizedPower));
end

