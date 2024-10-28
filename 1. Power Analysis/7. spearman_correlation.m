function run_spearman_correlation_analysis()
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % List of relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

       controlAges = containers.Map({'71075', '71075a', '75872', '75872a', '77862', '77862a', '83810', ...
                                  '87925', '92359', '110903', '113566', '113566a', '119655', '121139', ...
                                  '121139a', '144955', '144955a', '158439', '158439a', '171610', ...
                                  '171679', '172404'}, ...
                                  [6.5, 6.5, 9.1, 9.1, 12.4, 12.4, 9.8, 7.5, 8.2, 9.5, 10, 10, 8, ...
                                  9.2, 9.2, 11.7, 11.7, 13.8, 13.8, 8.1, 7.8, 7.6]);
    
    dsAges = containers.Map({'23281', '73716', '79426', '92187', '92446', '111822', ...
                             '138825', '148074', '167741'}, ...
                             [15.11, 9.1, 9.1, 8.6, 10, 16.5, 6.1, 16.2, 11.3]);
    % Process files for the Control group
    fprintf('Processing Control EEG files...\n');
    [controlDeltaPower, controlThetaPower, agesControl] = process_files(controlFolderPath, relevantElectrodes, controlAges);

    % Process files for the DS group
    fprintf('Processing DS EEG files...\n');
    [dsDeltaPower, dsThetaPower, agesDS] = process_files(dsFolderPath, relevantElectrodes, dsAges);

    % Perform Spearman Correlation Analysis
    fprintf('Performing Spearman Correlation Analysis...\n');
    try
        % Control group delta and theta
        [rhoDeltaControl, pDeltaControl] = perform_spearman_correlation(controlDeltaPower, agesControl);
        [rhoThetaControl, pThetaControl] = perform_spearman_correlation(controlThetaPower, agesControl);
        fprintf('Control Group Delta Power: rho = %.4f, p = %.4f\n', rhoDeltaControl, pDeltaControl);
        fprintf('Control Group Theta Power: rho = %.4f, p = %.4f\n', rhoThetaControl, pThetaControl);

        % DS group delta and theta
        [rhoDeltaDS, pDeltaDS] = perform_spearman_correlation(dsDeltaPower, agesDS);
        [rhoThetaDS, pThetaDS] = perform_spearman_correlation(dsThetaPower, agesDS);
        fprintf('DS Group Delta Power: rho = %.4f, p = %.4f\n', rhoDeltaDS, pDeltaDS);
        fprintf('DS Group Theta Power: rho = %.4f, p = %.4f\n', rhoThetaDS, pThetaDS);
    catch ME
        fprintf('Error during Spearman Correlation Analysis: %s\n', ME.message);
    end
end

%% Helper Function: Process files and compute δ and θ power using Welch's method
function [deltaPower, thetaPower, ages] = process_files(folderPath, electrodes, ageMap)
    files = dir(fullfile(folderPath, '*.edf'));
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
                        if iscell(signal)
                            signal = cell2mat(signal);
                        end
                        electrodeData = [electrodeData, signal(~isnan(signal))];
                    else
                        fprintf('Warning: Electrode %s not found in file %s\n', electrode, files(k).name);
                    end
                end

                % Compute δ (delta) and θ (theta) power using Welch's method
                if ~isempty(electrodeData)
                    Fs = info.NumSamples(1); % Sampling frequency

                    % Calculate power spectrum with Welch's method
                    [frequencies, powerSpectrum] = compute_welch_power(electrodeData, Fs);

                    % Delta power (0.5 - 4 Hz)
                    deltaIndices = frequencies >= 0.5 & frequencies <= 4;
                    deltaPower = [deltaPower; mean(powerSpectrum(deltaIndices))];

                    % Theta power (4 - 8 Hz)
                    thetaIndices = frequencies >= 4 & frequencies <= 8;
                    thetaPower = [thetaPower; mean(powerSpectrum(thetaIndices))];

                    % Store age
                    age = ageMap(fileName);
                    ages = [ages; age];

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

    fprintf('Processed %d files: Delta Power Size: %d, Theta Power Size: %d, Ages Size: %d\n', length(files), length(deltaPower), length(thetaPower), length(ages));
end

%% Helper Function: Compute Welch Power Spectrum
function [frequencies, powerSpectrum] = compute_welch_power(signal, Fs)
    windowLength = 2 * Fs; % 2-second window
    overlap = 0.5 * windowLength; % 50% overlap
    [powerSpectrum, frequencies] = pwelch(signal, hann(windowLength), overlap, [], Fs);
end

%% Helper Function: Perform Spearman Correlation
function [rho, pValue] = perform_spearman_correlation(data1, data2)
    [rho, pValue] = corr(data1, data2, 'Type', 'Spearman');
end
