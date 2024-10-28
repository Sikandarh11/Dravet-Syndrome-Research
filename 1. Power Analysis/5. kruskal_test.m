function kruskal_test
    % Load precomputed power spectra from Welch's method
    load('power_spectra.mat', 'controlPower', 'dsPower');

    % Ensure control and DS power data are not empty
    if isempty(controlPower) || isempty(dsPower)
        error('Power spectra data not found. Ensure that Welch''s method has been executed and saved.');
    end

    % Ensure controlPower and dsPower are column vectors for comparison
    controlPower = reshape(controlPower, [], 1);
    dsPower = reshape(dsPower, [], 1);

    % Perform Kruskal-Wallis Test using Welch's power spectra data
    fprintf('Performing Kruskal-Wallis Test on Welch power data...\n');
    perform_kruskal_test(controlPower, dsPower);
end

%% Helper Function: Perform Kruskal-Wallis Test
function perform_kruskal_test(controlData, dsData)
    fprintf('Performing Kruskal-Wallis Test...\n');
    try
        % Combine data and group labels
        data = [controlData; dsData];
        groupLabels = [ones(length(controlData), 1); 2 * ones(length(dsData), 1)]; % 1 for control, 2 for DS

        % Perform Kruskal-Wallis test
        [p, tbl, stats] = kruskalwallis(data, groupLabels, 'off');
        
        % Display results
        fprintf('Kruskal-Wallis Test p-value: %.4f\n', p);
        disp(tbl);
        disp(stats);
    catch ME
        fprintf('Error in Kruskal-Wallis Test: %s\n', ME.message);
    end
end
