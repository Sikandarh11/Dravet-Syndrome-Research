function pac_subanalysis()
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
    
    % Segment duration parameters
    numSegments = 2;  % First divide data into halves
    numThirds = 3;     % Then divide data into thirds

    % Debugging flag to check each step
    debug = true;

    % Process Control and DS EEG files for PAC
    fprintf('Processing Control EEG files for PAC subanalysis...\n');
    controlPAC = process_pac_files(controlFolderPath, controlAges, 0, numSegments, numThirds, debug);
    
    fprintf('Processing DS EEG files for PAC subanalysis...\n');
    dsPAC = process_pac_files(dsFolderPath, dsAges, 1, numSegments, numThirds, debug);

    % Concatenate Control and DS data for ANOVA
    fprintf('Running ANOVA on PAC results...\n');
    allPACData = [controlPAC; dsPAC];
    
    % ANOVA on segments and group differences
    run_anova_and_permutation_tests(allPACData, numSegments, numThirds);
end

%% Helper Function: Process Files and Perform PAC Subanalysis
function pacData = process_pac_files(folderPath, ageMap, label, numSegments, numThirds, debug)
    files = dir(fullfile(folderPath, '*.edf'));
    pacData = [];

    for k = 1:length(files)
        try
            edfFile = fullfile(folderPath, files(k).name);
            data = edfread(edfFile);

            % Extract data for relevant electrodes (assuming predefined electrodes)
            signalData = extract_electrode_data(data);
            
            % Split into halves and thirds
            [halfSegments, thirdSegments] = segment_data(signalData, numSegments, numThirds, debug);
            
            % Calculate PAC for each segment
            pacResults = calculate_pac_for_segments(halfSegments, thirdSegments, edfFile, debug);

            % Extract age from file name
            fileID = extractBefore(files(k).name, '.edf');
            if isKey(ageMap, fileID)
                age = ageMap(fileID);
                for resultIdx = 1:length(pacResults)
                    pacData = [pacData; pacResults(resultIdx), age, label, resultIdx]; %#ok<AGROW>
                end
            else
                fprintf('Age not found for file: %s\n', files(k).name);
            end
        catch ME
            fprintf('Error processing %s: %s\n', files(k).name, ME.message);
        end
    end
end

%% Helper Function: Segment Data into Halves and Thirds
function [halfSegments, thirdSegments] = segment_data(signalData, numSegments, numThirds, debug)
    % Segment the data into halves and thirds
    segmentLengthHalf = floor(length(signalData) / numSegments);
    segmentLengthThird = floor(length(signalData) / numThirds);

    % Divide the signal into halves
    halfSegments = cell(numSegments, 1);
    for i = 1:numSegments
        halfSegments{i} = signalData((i-1)*segmentLengthHalf + 1 : i*segmentLengthHalf);
        if debug
            fprintf('Half Segment %d: Length = %d\n', i, length(halfSegments{i}));
        end
    end

    % Divide the signal into thirds
    thirdSegments = cell(numThirds, 1);
    for i = 1:numThirds
        thirdSegments{i} = signalData((i-1)*segmentLengthThird + 1 : i*segmentLengthThird);
        if debug
            fprintf('Third Segment %d: Length = %d\n', i, length(thirdSegments{i}));
        end
    end
end

%% Helper Function: Calculate PAC for each segment
function pacResults = calculate_pac_for_segments(halfSegments, thirdSegments, fileName, debug)
    pacResults = [];

    % Calculate PAC for halves
    for i = 1:length(halfSegments)
        pacMI = calculate_pac(halfSegments{i});
        pacResults = [pacResults; pacMI]; %#ok<AGROW>
        if debug
            fprintf('Processed Half Segment %d of %s: PAC MI = %.4f\n', i, fileName, pacMI);
        end
    end

    % Calculate PAC for thirds
    for i = 1:length(thirdSegments)
        pacMI = calculate_pac(thirdSegments{i});
        pacResults = [pacResults; pacMI]; %#ok<AGROW>
        if debug
            fprintf('Processed Third Segment %d of %s: PAC MI = %.4f\n', i, fileName, pacMI);
        end
    end
end

%% Helper Function: Run ANOVA and Permutation Tests
function run_anova_and_permutation_tests(pacData, numSegments, numThirds)
    % Prepare data for ANOVA: Group (Control/DS) and time segments (halves/thirds)
    groupLabels = pacData(:, 3);  % Control = 0, DS = 1
    ages = pacData(:, 2);         % Participant ages
    pacMI = pacData(:, 1);        % PAC values

    % ANOVA: Investigating group differences and time segment stability
    fprintf('Running Two-Way Repeated Measures ANOVA...\n');
    groupFactor = groupLabels;
    timeFactor = pacData(:, 4);  % Time segment index (half/third)

    % Perform ANOVA on the data
    [p, tbl, stats] = anovan(pacMI, {groupFactor, timeFactor}, 'model', 'interaction', ...
                             'varnames', {'Group', 'TimeSegment'});

    % Print ANOVA results
    fprintf('ANOVA Results:\n');
    disp(tbl);
    
    % Posthoc: Multiple comparisons and permutation testing
    fprintf('Running Posthoc Multiple Comparisons and Permutation Tests...\n');
    perform_posthoc_tests(pacData, groupLabels, numSegments, numThirds);
end

%% Helper Function: Posthoc Multiple Comparisons and Permutation Testing
function perform_posthoc_tests(pacData, groupLabels, numSegments, numThirds)
    % Split data into halves and thirds for multiple comparisons
    fprintf('Posthoc Multiple Comparisons and Permutation Tests...\n');
    
    % Permutation testing (randomly shuffle group labels and compare)
    numPermutations = 1000;
    pValues = [];
    for i = 1:numPermutations
        shuffledLabels = groupLabels(randperm(length(groupLabels)));
        % Perform independent t-test for shuffled groups
        [~, pValue] = ttest2(pacData(groupLabels == 0, 1), pacData(groupLabels == 1, 1));
        pValues = [pValues; pValue];
        fprintf('Permutation %d: p-value = %.4f\n', i, pValue);
    end
    
    % Apply Bonferroni correction for multiple comparisons
    correctedPValues = min(pValues * numPermutations, 1);
    fprintf('Bonferroni Corrected p-values:\n');
    disp(correctedPValues);
end

%% Placeholder for PAC Calculation (to be modified based on actual method)
function pacMI = calculate_pac(segment)
    % Perform PAC calculation for a given segment using the modulation index (MI) method
    % This is a placeholder for the actual PAC calculation (e.g., Tort MI, KL MI)
    pacMI = abs(mean(hilbert(segment)));  % Using Hilbert as placeholder
end

%% Helper Function: Extract Data for Relevant Electrodes
function electrodeData = extract_electrode_data(data)
    % List of relevant electrodes (modify this based on the actual channels in your data)
    relevantElectrodes = ["EEGF3_Cz", "EEGFz_Cz", "EEGF4_Cz", ...
                          "EEGC3_Cz", "EEGCz_Cz", "EEGC4_Cz", ...
                          "EEGP3_Cz", "EEGPz_Cz", "EEGP4_Cz"];

    % Initialize storage for electrode data
    electrodeData = [];

    % Extract relevant electrode data from timetable
    for electrode = relevantElectrodes
        if ismember(electrode, data.Properties.VariableNames)
            % Extract signal data for this electrode
            signal = data{:, electrode};

            % Ensure the data is numeric and remove NaN values
            signal = cell2mat(signal); % Convert cell array to numeric array if needed
            signal = signal(~isnan(signal)); % Remove NaN values

            % Append the cleaned signal to the electrodeData matrix
            electrodeData = [electrodeData, signal]; %#ok<AGROW>
        else
            fprintf('Electrode %s not found in file.\n', electrode);
        end
    end

    % Check if any electrode data was found, otherwise return empty
    if isempty(electrodeData)
        error('No valid electrode data found for the specified electrodes.');
    end
end
