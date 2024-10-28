function pac_analysis
    % Define folder paths for Control EEG and DS cases
    controlFolderPath = '/MATLAB Drive/EEG newDataset/Control EEG/';
    dsFolderPath = '/MATLAB Drive/EEG newDataset/DS cases/';

    % Control and DS age maps
    controlAges = containers.Map({'71075', '71075a', '75872', '75872a', '77862', '77862a', '83810', ...
                                  '87925', '92359', '110903', '113566', '113566a', '119655', '121139', ...
                                  '121139a', '144955', '144955a', '158439', '158439a', '171610', ...
                                  '171679', '172404'}, ...
                                  [6.5, 6.5, 9.1, 9.1, 12.4, 12.4, 9.8, 7.5, 8.2, 9.5, 10, 10, 8, ...
                                  9.2, 9.2, 11.7, 11.7, 13.8, 13.8, 8.1, 7.8, 7.6]);

    dsAges = containers.Map({'23281', '73716', '79426', '92187', '92446', '111822', ...
                             '138825', '148074', '167741'}, ...
                             [15.11, 9.1, 9.1, 8.6, 10, 16.5, 6.1, 16.2, 11.3]);

    % Relevant electrodes
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % PAC calculation parameters
    numBins = 18;  % Phase bins for PAC calculation

    % Process Control and DS EEG files
    fprintf('Processing Control EEG files...\n');
    controlData = process_files_for_pac(controlFolderPath, controlAges, relevantElectrodes, numBins, 0);

    fprintf('Processing DS EEG files...\n');
    dsData = process_files_for_pac(dsFolderPath, dsAges, relevantElectrodes, numBins, 1);

    % Prepare data for logistic regression and statistical tests
    allData = [controlData.Data; dsData.Data];
    allLabels = [controlData.Labels; dsData.Labels];
    ages = [controlData.Ages; dsData.Ages];

    % Run Logistic Regression on PAC data
    fprintf('Running logistic regression on PAC data...\n');
    run_logistic_regression(allData, allLabels, ages);

    % Perform Statistical Analysis
    perform_statistical_tests(controlData.Data, dsData.Data, ages, controlData.Ages, dsData.Ages);
end

%% Helper Function: Process Files and Extract Data for PAC
function groupData = process_files_for_pac(folderPath, ageMap, relevantElectrodes, numBins, label)
    files = dir(fullfile(folderPath, '*.edf'));
    groupData.Data = [];
    groupData.Labels = [];
    groupData.Ages = [];

    for k = 1:length(files)
        try
            edfFile = fullfile(folderPath, files(k).name);
            data = edfread(edfFile);

            % Extract data for relevant electrodes
            electrodeData = [];
            for electrode = relevantElectrodes
                if ismember(electrode, data.Properties.VariableNames)
                    signal = data{:, electrode};
                    signal = cell2mat(signal); 
                    signal = signal(~isnan(signal));
                    electrodeData = [electrodeData, signal];
                end
            end

            % PAC calculation
            MI = calculate_pac_mi(electrodeData, numBins);

            % Extract age from file name
            fileID = extractBefore(files(k).name, '.edf');
            if isKey(ageMap, fileID)
                age = ageMap(fileID);
                groupData.Ages = [groupData.Ages; age];
                groupData.Data = [groupData.Data; MI];
                groupData.Labels = [groupData.Labels; label];
                fprintf('Processed %s with PAC MI: %.4f, Age: %.2f\n', files(k).name, MI, age);
            else
                fprintf('Age not found for file: %s\n', files(k).name);
            end
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end
end

%% Helper Function: PAC Calculation Using MI (Modulation Index)
function MI = calculate_pac_mi(electrodeData, numBins)
    analyticSignal = hilbert(electrodeData);
    phase = angle(analyticSignal);
    amplitude = abs(analyticSignal);

    % Bin phase data
    binEdges = linspace(-pi, pi, numBins + 1);
    meanAmplitude = zeros(1, numBins);

    % Mean amplitude in each bin
    for binIdx = 1:numBins
        binIndices = (phase >= binEdges(binIdx)) & (phase < binEdges(binIdx + 1));
        meanAmplitude(binIdx) = mean(amplitude(binIndices), 'omitnan');
    end

    % Normalize amplitudes and compute MI
    normalizedAmplitude = meanAmplitude / sum(meanAmplitude);
    uniformDist = 1 / numBins;
    log_N = log(numBins);
    MI = (log_N + sum(normalizedAmplitude .* log(normalizedAmplitude / uniformDist))) / log_N;
end

%% Helper Function: Logistic Regression
function run_logistic_regression(allData, allLabels, ages)
    predictorTable = table(allData, ages, allLabels, ...
        'VariableNames', {'PAC_MI', 'Age', 'GroupLabel'});

    fprintf('Fitting Logistic Regression Model...\n');
    glmModel = fitglm(predictorTable, 'GroupLabel ~ 1 + PAC_MI + Age', ...
                      'Distribution', 'binomial', 'Link', 'logit');
    disp(glmModel);

    % Prediction and accuracy
    fprintf('Predicting group membership...\n');
    predictedProbs = predict(glmModel, predictorTable);
    predictedLabels = predictedProbs >= 0.5;
    accuracy = mean(predictedLabels == allLabels) * 100;
    fprintf('Classification Accuracy: %.2f%%\n', accuracy);
end

%% Helper Function: Perform Statistical Tests
function perform_statistical_tests(controlData, dsData, allAges, controlAges, dsAges)
    % Custom Shapiro-Wilk Test for Normality (Alternative to `swtest`)
    fprintf('Running Shapiro-Wilk Test for Normality...\n');
    [hControl, pControl] = custom_shapiro_wilk(controlData);
    [hDS, pDS] = custom_shapiro_wilk(dsData);
    fprintf('Shapiro-Wilk Test - Control Group: p=%.4f, DS Group: p=%.4f\n', pControl, pDS);

    % Two-sample F-test for Equality of Variances
    [hVar, pVar] = vartest2(controlData, dsData);
    fprintf('Two-sample F-test for Variance Equality: p=%.4f\n', pVar);

    % Two-sample t-test for Mean Comparison (Bonferroni corrected)
    [hTtest, pTtest] = ttest2(controlData, dsData, 'Vartype', 'unequal');
    bonferroniCorrectedP = pTtest * 2; % Assuming only two comparisons
    fprintf('Two-sample t-test (unpaired): p=%.4f, Bonferroni corrected p=%.4f\n', pTtest, bonferroniCorrectedP);

    % Spearman Correlation with Age
    [rhoControl, pControlAge] = corr(controlAges, controlData, 'Type', 'Spearman');
    [rhoDS, pDSAge] = corr(dsAges, dsData, 'Type', 'Spearman');
    fprintf('Spearman Correlation - Control Group: rho=%.4f, p=%.4f\n', rhoControl, pControlAge);
    fprintf('Spearman Correlation - DS Group: rho=%.4f, p=%.4f\n', rhoDS, pDSAge);
end

%% Custom Shapiro-Wilk Test for Normality
function [h, p] = custom_shapiro_wilk(data)
    % Approximate Shapiro-Wilk test logic; consider using MATLAB statistical toolbox for full functionality.
    % Output h is 0 if normality is not rejected, 1 otherwise.
    % Output p is the p-value of the test.
    n = length(data);
    sortedData = sort(data);
    expected = norminv(((1:n)' - 0.375) / (n + 0.25));
    a = (expected' * expected) \ expected' * sortedData;
    w = sum((sortedData - a * expected).^2) / sum((sortedData - mean(sortedData)).^2);
    p = 1 - chi2cdf(w, n - 1);
    h = p < 0.05; % Reject normality at 5% significance level
end
