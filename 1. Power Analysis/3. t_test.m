function t_test
    % Load the Welch's method results (control and DS power data)
    load('power_spectra.mat', 'controlPower', 'dsPower');

    % Perform Independent-Samples T-Test
    fprintf('Performing Independent-Samples T-Test on Welch Power Spectra...\n');
    [h, pValue, ci, stats] = independent_ttest(controlPower, dsPower);

    % Display Results
    fprintf('\nT-Test Results (Control vs DS cases based on Welch Power Spectra):\n');
    fprintf('Hypothesis test result (h): %d\n', h);
    fprintf('p-value: %.4f\n', pValue);
    fprintf('95%% Confidence Interval: [%.4f, %.4f]\n', ci(1), ci(2));
    fprintf('T-statistic: %.4f\n', stats.tstat);
    fprintf('Degrees of Freedom: %.2f\n', stats.df);
end

%% Helper Function: Perform Independent-Samples T-Test
function [h, pValue, ci, stats] = independent_ttest(group1, group2)
    % Perform two-sample T-test with unequal variances
    [h, pValue, ci, stats] = ttest2(group1(:), group2(:), 'Vartype', 'unequal');
end

