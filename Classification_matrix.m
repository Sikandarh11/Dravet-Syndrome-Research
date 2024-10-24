function improved_logistic_regression()
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % control and ds age maps 
    controlAges = containers.Map({'71075', '71075a', '75872', '75872a', '77862', '77862a', '83810', ...
                                  '87925', '92359', '110903', '113566', '113566a', '119655', '121139', ...
                                  '121139a', '144955', '144955a', '158439', '158439a', '171610', ...
                                  '171679', '172404'}, ...
                                  [6.5, 6.5, 9.1, 9.1, 12.4, 12.4, 9.8, 7.5, 8.2, 9.5, 10, 10, 8, ...
                                  9.2, 9.2, 11.7, 11.7, 13.8, 13.8, 8.1, 7.8, 7.6]);

    dsAges = containers.Map({'23281', '73716', '79426', '92187', '92446', '111822', ...
                             '138825', '148074', '167741'}, ...
                             [15.11, 9.1, 9.1, 8.6, 10, 16.5, 6.1, 16.2, 11.3]);

    % Relevant electrodes from the paper
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Process Control and DS EEG files for delta and theta power and PAC
    fprintf('Processing Control EEG files...\n');
    controlData = process_files_for_power_and_pac(controlFolderPath, controlAges, relevantElectrodes, 0); % Label 0 for Control

    fprintf('Processing DS EEG files...\n');
    dsData = process_files_for_power_and_pac(dsFolderPath, dsAges, relevantElectrodes, 1); % Label 1 for DS



    fprintf('Running separate logistic regressions...\n');
    
    % Delta Power Logistic Regression
    run_individual_logistic_regression(controlData, dsData, 'DeltaPower');
    
    % Theta Power Logistic Regression
    run_individual_logistic_regression(controlData, dsData, 'ThetaPower');
    
    % PAC Logistic Regression
    run_individual_logistic_regression(controlData, dsData, 'PAC');

    % Generating plots
    plot_power_graphs(controlData, dsData);
end

%% Helper Function: Process Files and Extract Data for Delta Power, Theta Power, and PAC
function groupData = process_files_for_power_and_pac(folderPath, ageMap, relevantElectrodes, label)
    files = dir(fullfile(folderPath, '*.edf'));
    groupData.DeltaPower = [];
    groupData.ThetaPower = [];
    groupData.PAC = [];
    groupData.Labels = [];
    groupData.Ages = [];

    for k = 1:length(files)
        try
            edfFile = fullfile(folderPath, files(k).name);
            data = edfread(edfFile);

            % Extract data for relevant electrodes
            electrodeData = [];
            minLength = inf; % Initialize with a large value to find the minimum length
            for electrode = relevantElectrodes
                if ismember(electrode, data.Properties.VariableNames)
                    signal = data{:, electrode};
                    signal = cell2mat(signal); % Convert to numeric array
                    signal = signal(~isnan(signal)); % Remove NaNs

                    % Update minLength to ensure all signals have the same length
                    if length(signal) < minLength
                        minLength = length(signal);
                    end

                    electrodeData = [electrodeData, signal]; %#ok<AGROW>
                end
            end

            % Check if valid electrode data is found
            if isempty(electrodeData)
                fprintf('No valid electrode data found in file: %s\n', files(k).name);
                continue;
            end

            % Truncate all electrode signals to the minimum length
            electrodeData = electrodeData(1:minLength, :);


            deltaPower = calculate_band_power(electrodeData, [1 4], 250); % δ (1-4 Hz)
            thetaPower = calculate_band_power(electrodeData, [4 7], 250); % θ (4-7 Hz)

            % Store the mean δ and θ power and PAC
            pacValue = calculate_pac(electrodeData); % Placeholder for PAC calculation
            groupData.DeltaPower = [groupData.DeltaPower; mean(deltaPower)];
            groupData.ThetaPower = [groupData.ThetaPower; mean(thetaPower)];
            groupData.PAC = [groupData.PAC; pacValue];

            % Extract age from the file name
            fileID = extractBefore(files(k).name, '.edf'); % Extract file ID
            if isKey(ageMap, fileID)
                age = ageMap(fileID);
                groupData.Ages = [groupData.Ages; age];
                groupData.Labels = [groupData.Labels; label]; % 0 for Control, 1 for DS
                fprintf('Processed %s with %d relevant signals. Age: %.2f, δ Power: %.2f, θ Power: %.2f, PAC: %.2f\n', ...
                        files(k).name, size(electrodeData, 2), age, mean(deltaPower), mean(thetaPower), pacValue);
            else
                fprintf('Age not found for file: %s\n', files(k).name);
            end
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end
end



%% Helper Function: Calculate Band Power
%%==================================== Question : Mam Soyiba===============
% for calculating band power i took the mean of the filtered data, bcz
% please review this part
%I used the process with or without mean, results are same
%I uploaded the results screenshots please check the results also
function bandPower = calculate_band_power(electrodeData, bandRange, Fs)
    % Apply band-pass filter for the specific frequency range
    [b, a] = butter(4, bandRange / (Fs / 2), 'bandpass');
    filteredData = filtfilt(b, a, electrodeData);
    disp('band  Data Values:');
    
    bandPower = mean(filteredData.^2, 1); % Power is the mean squared signal =======================================Question=============================================
    disp(bandpower);  
end

%% Helper Function: Calculate Phase-Amplitude Coupling (PAC)
function pacValue = calculate_pac(electrodeData)
    % Placeholder for real PAC calculation logic based on paper's method
    pacValue = abs(mean(hilbert(electrodeData(:)))); %%=======================================================Also please take a look at here=============================
end

%% Function to Run Logistic Regression for Specific Metric (Delta, Theta, or PAC)
function run_individual_logistic_regression(controlData, dsData, metric)

    % Prepare for logistic regression
    fprintf('Running logistic regression for %s...\n', metric);
    allData = [controlData.(metric); dsData.(metric)];
    allLabels = [controlData.Labels; dsData.Labels];
    ages = [controlData.Ages; dsData.Ages];

    
    % Z-Score Normalization
    fprintf('Applying Z-score normalization for %s and Age...\n', metric);
    
    % Normalize the metric (DeltaPower, ThetaPower, PAC)
    meanMetric = mean(allData);
    stdMetric = std(allData);
    normalizedMetric = (allData - meanMetric) / stdMetric;
    
    % Normalize Age
    meanAge = mean(ages);
    stdAge = std(ages);
    normalizedAge = (ages - meanAge) / stdAge;

    % Prepare predictor table for GLM with normalized values
    predictorTable = table(normalizedMetric, normalizedAge, allLabels, 'VariableNames', {metric, 'Age', 'GroupLabel'});

    % Fittinh GLM model with a binomial distribution
    try
        glmModel = fitglm(predictorTable, sprintf('GroupLabel ~ 1 + %s + Age', metric), ...
                          'Distribution', 'binomial', 'Link', 'logit');
        disp(glmModel);
    catch ME
        fprintf('Error fitting GLM model for %s: %s\n', metric, ME.message);
        return;
    end

    % Predict probabilities for each participant
    fprintf('Predicting group membership for %s...\n', metric);
    predictedProbs = predict(glmModel, predictorTable);

    % Compute accuracy
    predictedLabels = predictedProbs >= 0.5; % Threshold at 0.5
    correctPredictions = predictedLabels == allLabels;
    accuracy = mean(correctPredictions) * 100;

    % Display results
    fprintf('%s Classification Accuracy: %.2f%%\n', metric, accuracy);
end




%% Helper Function: Plot Power Graphs
function plot_power_graphs(controlData, dsData)
    % Create plots for δ and θ power distributions and logistic regression results
    
    % Plot relative power for delta (1-4 Hz) and theta (4-7 Hz)
    figure;
    
    % Plot Delta Power Distribution (Figure A)
    subplot(2, 2, 1); 
    hold on;
    plot_relative_power_distribution(controlData.DeltaPower, 'Control');
    plot_relative_power_distribution(dsData.DeltaPower, 'DS');
    title('Relative Delta Power (1-4 Hz)');
    legend('Control', 'DS');
    xlabel('Power');
    ylabel('Density');
    hold off;
    
    % Plot Theta Power Distribution (Figure B)
    subplot(2, 2, 2);
    hold on;
    plot_relative_power_distribution(controlData.ThetaPower, 'Control');
    plot_relative_power_distribution(dsData.ThetaPower, 'DS');
    title('Relative Theta Power (4-7 Hz)');
    legend('Control', 'DS');
    xlabel('Power');
    ylabel('Density');
    hold off;

    % Boxplots for δ and θ power (Figures C and D)
    subplot(2, 2, 3);
    boxplot([controlData.DeltaPower; dsData.DeltaPower], [controlData.Labels; dsData.Labels], 'Labels', {'Control', 'DS'});
    title('Delta Power Comparison');
    ylabel('Power');
    
    subplot(2, 2, 4);
    boxplot([controlData.ThetaPower; dsData.ThetaPower], [controlData.Labels; dsData.Labels], 'Labels', {'Control', 'DS'});
    title('Theta Power Comparison');
    ylabel('Power');

    % Creating new figure for logistic regression scatter plots (Figures E and F)
    figure;
    subplot(1, 2, 1); % Delta logistic regression scatter plot
    scatter_logistic_results(controlData.DeltaPower, dsData.DeltaPower, 'Delta Power');
    
    subplot(1, 2, 2); % Theta logistic regression scatter plot
    scatter_logistic_results(controlData.ThetaPower, dsData.ThetaPower, 'Theta Power');
end

%% Helper Function: Plot Relative Power Distribution
function plot_relative_power_distribution(powerData, label)
    % Helper to plot power distribution curves
    [f, xi] = ksdensity(powerData); % Kernel density estimation for smooth curve
    plot(xi, f, 'DisplayName', label);
end

%% Helper Function: Scatter Logistic Regression Results
function scatter_logistic_results(controlData, dsData, titleText)
    % Create scatter plot for logistic regression results (Figures E and F)
    scatter(ones(size(controlData)), controlData, 'r', 'filled');
    hold on;
    scatter(2*ones(size(dsData)), dsData, 'b', 'filled');
    title([titleText ' Classification']);
    ylabel('Power');
    xticks([1 2]);
    xticklabels({'Control', 'DS'});
    legend('Control', 'DS');
    hold off;
end
