function bonferroni_correction
    % Load precomputed power spectra from Welch's method
    load('power_spectra.mat', 'controlPower', 'dsPower');

    % Define frequency ranges: delta (1-4 Hz), theta (4-8 Hz), alpha (8-13 Hz), beta (13-30 Hz)
    frequencyRanges = [1 4; 4 8; 8 13; 13 30]; 

    % Initialize storage for p-values
    pValues = []; 

    % Perform T-Test for each frequency range
    fprintf('Comparing frequency ranges with T-Test...\n');
    for i = 1:size(frequencyRanges, 1)
        fprintf('Testing frequency range: %d-%d Hz\n', frequencyRanges(i, 1), frequencyRanges(i, 2));

        % Extract frequency-specific data from control and DS power data
        controlFreqData = extract_frequency_data(controlPower, frequencyRanges(i, :));
        dsFreqData = extract_frequency_data(dsPower, frequencyRanges(i, :));

        % Perform T-Test
        [~, pValue] = ttest2(controlFreqData, dsFreqData, 'Vartype', 'unequal');
        pValues = [pValues; pValue];

        fprintf('p-value for range %d-%d Hz: %.4f\n', frequencyRanges(i, 1), frequencyRanges(i, 2), pValue);
    end

    % Apply Bonferroni correction
    fprintf('Applying Bonferroni Correction...\n');
    alpha = 0.05;
    [adjustedPValues, significant] = apply_bonferroni_correction(pValues, alpha);

    % Display results
    fprintf('Adjusted p-values:\n');
    disp(adjustedPValues);
    fprintf('Significance Status (1 = significant, 0 = not significant):\n');
    disp(significant);
end

%% Helper Function: Apply Bonferroni Correction
function [adjustedPValues, significant] = apply_bonferroni_correction(pValues, alpha)
    % Adjusts p-values for multiple comparisons using Bonferroni correction
    m = length(pValues); 
    adjustedPValues = min(pValues * m, 1);
    significant = adjustedPValues < alpha;
end

%% Helper Function: Extract Frequency Data Based on Range
function freqData = extract_frequency_data(powerData, freqRange)
    % Ensure we stay within the data bounds
    freqRange = min(freqRange(1):min(freqRange(2), size(powerData, 2)));
    % Average power within the range for each entry
    freqData = mean(powerData(:, freqRange), 2); 
end
