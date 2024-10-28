function f_test_on_welch
    % Load the Welch's method results
    load('power_spectra.mat', 'controlPower', 'dsPower');

    % Perform F-Test on the normalized power data
    fprintf('Performing F-Test on Welch Power Spectra...\n');
    [h, pValue, F, varControl, varDS] = ftest2(controlPower, dsPower);

    % Display Results
    fprintf('F-Test Results:\n');
    fprintf('Hypothesis test result (h): %d\n', h);
    fprintf('p-value: %.4f\n', pValue);
    fprintf('F-statistic: %.4f\n', F);
    fprintf('Control group variance: %.4f\n', varControl);
    fprintf('DS group variance: %.4f\n', varDS);
end

%% Helper Function: F-Test Calculation
function [h, pValue, F, var1, var2] = ftest2(group1, group2)
    % Calculate variances
    var1 = var(group1(:));
    var2 = var(group2(:));

    % Calculate F-statistic
    F = var1 / var2;
    n1 = length(group1) - 1;
    n2 = length(group2) - 1;

    % Two-tailed p-value
    pValue = 2 * (1 - fcdf(max(F, 1/F), n1, n2));

    % Hypothesis result: true if p-value < 0.05
    h = pValue < 0.05;
end
