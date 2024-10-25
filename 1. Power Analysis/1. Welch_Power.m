function welch_power
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % File patterns to match EDF files
    controlFilePattern = fullfile(controlFolderPath, '*.edf');
    dsFilePattern = fullfile(dsFolderPath, '*.edf');

    % Get the list of all EDF files in both directories
    controlFiles = dir(controlFilePattern);
    dsFiles = dir(dsFilePattern);

    % Debugging: Verify directories and files found
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

    % Output final results
    fprintf('Control Group Power Data:\n');
    disp(controlPower);

    fprintf('DS Group Power Data:\n');
    disp(dsPower);
end

%% Helper Function: Process Group Files and Extract Power Spectra
function powerData = process_group_files(fileList, folderPath)
    powerData = [];
    for k = 1:length(fileList)
        try
            filePath = fullfile(folderPath, fileList(k).name);
            info = edfinfo(filePath);
            data = edfread(filePath);

            % Extract EEG signal data
            signalData = data{:, 1};
            numericSignalData = cell2mat(signalData);
            signalData = numericSignalData(~isnan(numericSignalData));

            % Compute power spectrum using Welch's method
            [frequencies, power] = compute_welch_power(signalData);

            % Normalize power within the 1-80 Hz range
            normalizedPower = normalize_power(frequencies, power);
            powerData = [powerData; normalizedPower]; %#ok<AGROW>

            % Debugging: Output power spectrum status
            fprintf('Processed %s. Normalized Power: %s\n', ...
                fileList(k).name, mat2str(normalizedPower));
        catch ME
            fprintf('Error processing %s: %s\n', fileList(k).name, ME.message);
        end
    end
end

%% Helper Function: Compute Power Spectrum Using Welch’s Method
function [frequencies, power] = compute_welch_power(signal)
    % Parameters for Welch's method
    windowLength = 2 * 256;  % 2-second Hann window (assuming 256 Hz sample rate)
    overlap = 0.5 * windowLength;  % 50% overlap

    % Apply Welch's method to compute power spectrum
    [power, frequencies] = pwelch(signal, hann(windowLength), overlap, [], 256);

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



