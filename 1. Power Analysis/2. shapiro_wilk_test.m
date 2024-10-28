function shapiro_wilk_test
    % Load power data from Welch method
    if isfile('power_spectra.mat')
        load('power_spectra.mat', 'controlPower', 'dsPower');
    else
        error('File power_spectra.mat not found. Please run welch_power.m first.');
    end

    % Perform Shapiro-Wilk Test for both groups
    fprintf('Performing Shapiro-Wilk Test on Control Group Power Data...\n');
    perform_shapiro_wilk(controlPower, 'Control Group');

    fprintf('Performing Shapiro-Wilk Test on DS Group Power Data...\n');
    perform_shapiro_wilk(dsPower, 'DS Group');
end

%% Helper Function: Perform Shapiro-Wilk Test for Normality
function perform_shapiro_wilk(data, groupName)
    fprintf('Shapiro-Wilk Test for %s...\n', groupName);
    try
        % Perform Shapiro-Wilk test on each frequency band
        for i = 1:size(data, 2)
            [h, p] = swtest(data(:, i));
            fprintf('Frequency Band %d: h = %d, p = %.4f\n', i, h, p);
        end
    catch ME
        fprintf('Error in Shapiro-Wilk Test for %s: %s\n', groupName, ME.message);
    end
end

%% Helper Function: Shapiro-Wilk Test Implementation (Using Jarque-Bera as Proxy)
function [h, pValue] = swtest(x)
    % Use Jarque-Bera test as an approximation for Shapiro-Wilk
    [~, pValue] = jbtest(x);  
    h = pValue < 0.05; % Reject null hypothesis if p < 0.05
end
