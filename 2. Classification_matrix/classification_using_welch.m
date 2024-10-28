function classification()
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Process Control and DS EEG files using Welch's method to calculate power and PAC
    fprintf('Processing Control EEG files with Welch method...\n');
    controlData = process_files_with_welch(controlFolderPath, relevantElectrodes, 0); % Label 0 for Control

    fprintf('Processing DS EEG files with Welch method...\n');
    dsData = process_files_with_welch(dsFolderPath, relevantElectrodes, 1); % Label 1 for DS
    
    % Run logistic regression for each metric
    fprintf('Running logistic regressions for Delta Power, Theta Power, and PAC...\n');
    run_individual_logistic_regression(controlData, dsData, 'DeltaPower');
    run_individual_logistic_regression(controlData, dsData, 'ThetaPower');
    run_individual_logistic_regression(controlData, dsData, 'PAC');
end

%% Helper Function: Process Files with Welch's Method for Delta Power, Theta Power, and PAC
function groupData = process_files_with_welch(folderPath, relevantElectrodes, label)
    files = dir(fullfile(folderPath, '*.edf'));
    groupData.DeltaPower = [];
    groupData.ThetaPower = [];
    groupData.PAC = [];
    groupData.Labels = [];

    for k = 1:length(files)
        try
            edfFile = fullfile(folderPath, files(k).name);
            data = edfread(edfFile);

            % Extract data for relevant electrodes
            electrodeData = [];
            for electrode = relevantElectrodes
                if ismember(electrode, data.Properties.VariableNames)
                    signal = cell2mat(data{:, electrode});
                    electrodeData = [electrodeData, signal(~isnan(signal))];
                end
            end

            % Ensure we have data for analysis
            if isempty(electrodeData)
                fprintf('No valid electrode data found in file: %s\n', files(k).name);
                continue;
            end

            % Compute Welch power spectrum
            deltaPower = calculate_welch_power(electrodeData, [1 4], 250); % δ (1-4 Hz)
            thetaPower = calculate_welch_power(electrodeData, [4 7], 250); % θ (4-7 Hz)
            pacValue = calculate_pac(electrodeData); % Placeholder for PAC calculation
            
            % Store the results
            groupData.DeltaPower = [groupData.DeltaPower; mean(deltaPower)];
            groupData.ThetaPower = [groupData.ThetaPower; mean(thetaPower)];
            groupData.PAC = [groupData.PAC; pacValue];
            groupData.Labels = [groupData.Labels; label];

            fprintf('Processed %s with Delta Power: %.2f, Theta Power: %.2f, PAC: %.2f\n', ...
                files(k).name, mean(deltaPower), mean(thetaPower), pacValue);
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end
end

%% Helper Function: Calculate Welch Power Spectrum
function power = calculate_welch_power(electrodeData, bandRange, Fs)
    % Apply band-pass filter for the specific frequency range
    [b, a] = butter(4, bandRange / (Fs / 2), 'bandpass');
    filteredData = filtfilt(b, a, electrodeData);
    power = mean(filteredData .^ 2); % Power as mean squared signal
end

%% Helper Function: Placeholder for PAC Calculation
function pacValue = calculate_pac(electrodeData)
    % Placeholder for real PAC calculation logic based on paper's method
    pacValue = abs(mean(hilbert(electrodeData(:))));
end

%% Function to Run Logistic Regression for Each Metric
function run_individual_logistic_regression(controlData, dsData, metric)
    % Prepare data for logistic regression
    fprintf('Running logistic regression for %s...\n', metric);
    allData = [controlData.(metric); dsData.(metric)];
    allLabels = [controlData.Labels; dsData.Labels];

    % Z-Score Normalization
    meanMetric = mean(allData);
    stdMetric = std(allData);
    normalizedMetric = (allData - meanMetric) / stdMetric;

    % Prepare predictor table for GLM with normalized values
    predictorTable = table(normalizedMetric, allLabels, 'VariableNames', {metric, 'GroupLabel'});

    % Fit the GLM model with a binomial distribution
    glmModel = fitglm(predictorTable, sprintf('GroupLabel ~ 1 + %s', metric), ...
                      'Distribution', 'binomial', 'Link', 'logit');
    disp(glmModel);

    % Predict probabilities and compute accuracy
    predictedProbs = predict(glmModel, predictorTable);
    predictedLabels = predictedProbs >= 0.5; % Threshold at 0.5
    correctPredictions = predictedLabels == allLabels;
    accuracy = mean(correctPredictions) * 100;

    % Display results
    fprintf('%s Classification Accuracy: %.2f%%\n', metric, accuracy);
end
