%Tagtools substitute applicaiton 

%This uses the list of whale names created in previous scripts, if not yet
%applied, uncomment the first section 

%Uses the jay_find_dives function which will need to be on the workflow and
%the CircStats package which can be downloaded from https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
%And unzipped where Matlab can get to it

%The function SHOULD plot dive phases and append the dive shape and phase
%to each row. Do with this information as you will. 

%It also creates one megadataframe (c) with all the deployments in one
%object for easier holistic analysis

% Part 1: Load and Store whale_data Structures
%whale_names = dir('*processed.mat');

% Loop through the list and filter out files with underscores in their names
%for i = 1:length(whale_names)
    % Check if the file name does not contain an underscore
 %   if isempty(strfind(whale_names(i).name, '_'))
        % Load the file if it meets the criteria
%        data = load(whale_names(i).name);
        
        % Create a variable name without the .mat extension
 %       var_name = erase(whale_names(i).name, '.mat');
        
        % Assign the whale_data structure to the base workspace
%        assignin('base', var_name, data.whale_data);
%    end
%end
%%
%Change the directory instruction based on the naming
%conventions/classification tool you used
whale_names = dir("*processed_SVM_Classified.mat");
%Loop through the list and filter out files with underscores in their names
for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the whale name
    whale_name = whale_name(1:end-19); %20 for LSTM, 19 for SVM, if you have a different number of deployments change this number
    whale_data = load(whale_names(i).name); % Access the whale data object dynamicall
    % Save the trimmed data back into the whale object
    assignin('base', whale_name, whale_data.whale_data); % Update the object in the workspace
end
clear("i", "whale_name", "whale_data")
% Get all objects in the base workspace
whale_names = evalin('base', 'whos');

disp("Section 1, Data Loaded");
%%

tagtools_data = table();

for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    whale_data = evalin('base', whale_name);
    whale_data = jay_find_dives(whale_data, 5, 10, 1);
    saveas(gcf, ['Dive_Phases_', whale_name, '.png']);
    assignin('base', whale_name, whale_data); % Update the object in the workspace
end


% Initialize an empty table to hold data from all whale names

% Loop through each whale_name object and concatenate the specified columns
for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    current_data = evalin('base', whale_name);
    current_data.WhaleName = repmat(string(whale_name), height(current_data), 1);
    tagtools_data = [tagtools_data; current_data]; % Append to all_data table
end
close all; 

%%
%Summary Stats
% Filter the rows where MouthOpen == 1
mouth_open_data = tagtools_data(tagtools_data.MouthOpen == 1, :);

% Further filter by shape and class conditions
filtered_data = mouth_open_data(ismember(mouth_open_data.shape, ["U", "square"]) & mouth_open_data.class == "dest", :);

% Get the number of rows in the filtered data
total_filtered_rows = height(filtered_data);

% Get the total number of rows where MouthOpen == 1
total_mouth_open_rows = height(mouth_open_data);

% Calculate the percentage
percentage = (total_filtered_rows / total_mouth_open_rows) * 100;

% Display the result
fprintf('Percentage of Mouth Open time classified by diveshape: %.2f%%\n', percentage);

% Filter the rows satisfying the condition
filtered_data = tagtools_data(ismember(tagtools_data.shape, ["U", "square"]) & tagtools_data.class == "dest", :);

% Get the number of rows where MouthOpen == 1 in the filtered data
mouth_open_count = sum(filtered_data.MouthOpen == 1);

% Get the total number of rows in the filtered data
total_filtered_rows = height(filtered_data);

% Calculate the percentage
percentage = (mouth_open_count / total_filtered_rows) * 100;

% Display the result
fprintf('Percentage of dive shape rows identied as feeding with MouthOpen predicted: %.2f%%\n', percentage);

diveshape_count = height(tagtools_data(ismember(tagtools_data.shape, ["U", "square"]) & tagtools_data.class == "dest", :));
predicted_count = sum(tagtools_data.MouthOpen);
per = (predicted_count/diveshape_count)*100;
fprintf('Total dive shape rows as a percent of total predicted MouthOpen: %.2f%%\n', per);

%%
%Summary Stats by daylight classification
% Initialize a results table
results_table = table();

% Unique Daylight categories, including "Overall" for all data
daylight_categories = [unique(tagtools_data.Daylight); "Overall"];

% Loop over each Daylight category
for i = 1:length(daylight_categories)
    if daylight_categories(i) == "Overall"
        % Use the entire dataset for "Overall"
        current_data = tagtools_data;
    else
        % Filter the data for the current Daylight category
        current_data = tagtools_data(tagtools_data.Daylight == daylight_categories(i), :);
    end

    % Remove NaN values from the current dataset
    current_data = rmmissing(current_data);

    % Filter the rows where MouthOpen == 1
    mouth_open_data = current_data(current_data.MouthOpen == 1, :);

    % Further filter by shape and class conditions
    filtered_data = mouth_open_data(ismember(mouth_open_data.shape, ["U", "square"]) & mouth_open_data.class == "dest", :);

    % Calculate metrics for the current data
    total_filtered_rows = height(filtered_data);
    total_mouth_open_rows = height(mouth_open_data);
    percentage_mouth_open_dive_shape = (total_filtered_rows / total_mouth_open_rows) * 100;

    % Filter all rows satisfying the shape and class conditions
    filtered_data_all = current_data(ismember(current_data.shape, ["U", "square"]) & current_data.class == "dest", :);

    % Get the number of rows where MouthOpen == 1 in the filtered data
    mouth_open_count = sum(filtered_data_all.MouthOpen == 1);

    % Get the total number of rows in the filtered data
    total_filtered_rows_all = height(filtered_data_all);
    percentage_feeding_predicted = (mouth_open_count / total_filtered_rows_all) * 100;

    % Total rows with dive shape matching "U" or "square" and class "dest"
    diveshape_count = height(filtered_data_all);

    % Total predicted MouthOpen rows
    predicted_count = sum(current_data.MouthOpen);

    % Calculate percentage of dive shape rows identified as MouthOpen
    percentage_diveshape_mouth_open = (predicted_count / diveshape_count) * 100;

    % Append results for the current Daylight category to the table
    results_table = [results_table; table( ...
        daylight_categories(i), ...
        total_filtered_rows, ...
        total_mouth_open_rows, ...
        percentage_mouth_open_dive_shape, ...
        mouth_open_count, ...
        total_filtered_rows_all, ...
        percentage_feeding_predicted, ...
        diveshape_count, ...
        predicted_count, ...
        percentage_diveshape_mouth_open, ...
        'VariableNames', { ...
            'Daylight', ...
            'FilteredRowsMouthOpen', ...
            'TotalMouthOpenRows', ...
            'PercentMouthOpenDiveShape', ...
            'MouthOpenCount', ...
            'TotalFilteredRows', ...
            'PercentFeedingPredicted', ...
            'DiveShapeCount', ...
            'PredictedMouthOpenCount', ...
            'PercentDiveShapeMouthOpen'})];
end

% Display the results table
disp(results_table);
%%

%%
%Data combination
all_data = table(); % Initialize an empty table to hold data from all whale names

% Loop through each whale_name object and concatenate the specified columns
for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    current_data = evalin('base', whale_name);
    current_data.ID = repmat(string(whale_name), height(current_data), 1);
    all_data = [all_data; current_data]; % Append to all_data table
end
%%
%Mean daylight dest depth
mean_dest_day_depth = mean(all_data.Depth(strcmp(all_data.class, 'dest') & strcmp(all_data.Daylight, 'Day')));
disp(['Mean Depth for "dest" class during Day: ', num2str(mean_dest_day_depth)]);


std(all_data.Depth(strcmp(all_data.class, 'dest') & strcmp(all_data.Daylight, 'Day')));
%%
%Dive lengths 
% Initialize an empty table for dive lengths
dive_lengths = table();

for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    
    % Access the whale_data object dynamically
    whale_data = evalin('base', whale_name);
    
    % Get unique daylight values
    unique_daylight = unique(whale_data.Daylight);
    
    % Iterate through each daylight category
    for j = 1:length(unique_daylight)
        daylight_value = unique_daylight(j);
        
        % Filter data for current daylight category
        subset_data = whale_data(strcmp(whale_data.Daylight, daylight_value), :);
        
        % Find indices where 'class' is 'sur' within the subset
        sur_indices = find(strcmp(subset_data.class, 'sur'));
        
        % Identify non-consecutive occurrences
        gap_lengths = diff(sur_indices);  % Differences between consecutive 'sur' indices
        
        % Filter out cases where 'sur' is consecutive (gap == 1)
        non_consecutive_gaps = gap_lengths(gap_lengths > 1);
        
        % Compute the mean and standard deviation of gap lengths, handling empty cases
        if isempty(non_consecutive_gaps)
            mean_gap = 0; % Set to 0 if no non-consecutive gaps exist
            std_gap = 0;  % Standard deviation is also 0
        else
            mean_gap = mean(non_consecutive_gaps);
            std_gap = std(non_consecutive_gaps);
        end
        
        % Count the number of dives (new "sur" events)
        num_dives = sum([1; diff(sur_indices)] > 1); % Count starts of new dive events
        
        % Store results in the table, using a cell array for text variables
        % Divide the 10Hz data by 600 to get to minutes
        new_entry = table({whale_name}, {daylight_value}, mean_gap / 600, std_gap / 600, num_dives, ...
                          'VariableNames', {'deployment', 'Daylight', 'mean_gap', 'std_gap', 'num_dives'});
        
        % Append new row to the table
        dive_lengths = [dive_lengths; new_entry];
    end
end

dive_lengths.deployment = string(dive_lengths.deployment);
dive_lengths.Daylight = string(dive_lengths.Daylight);

% Display the results
disp(dive_lengths);

% Filter rows where Daylight is "Day"
day_indices = strcmp(dive_lengths.Daylight, "Day");

% Compute weighted mean of mean_gap for "Day" rows
mean_gap_day = sum(dive_lengths.mean_gap(day_indices) .* dive_lengths.num_dives(day_indices)) / sum(dive_lengths.num_dives(day_indices));

% Compute weighted standard deviation for "Day" rows
if sum(dive_lengths.num_dives(day_indices)) > 1
    std_gap_day = sqrt(sum(dive_lengths.num_dives(day_indices) .* (dive_lengths.std_gap(day_indices) .^ 2)) / sum(dive_lengths.num_dives(day_indices)));
else
    std_gap_day = 0; % If only one sample, standard deviation is 0
end

% Display the results
disp(['Mean gap during "Day": ', num2str(mean_gap_day)]);
disp(['Standard deviation of gap during "Day": ', num2str(std_gap_day)]);


%%
%Dive destination depths 
% Initialize arrays to store the results
day_mean_depths_all = [];
night_mean_depths_all = [];
day_std_depths_all = [];
night_std_depths_all = [];

% Group the data by deployment ID
deployment_ids = unique(all_data.ID);

for i = 1:length(deployment_ids)
    % Extract data for the current deployment
    current_deployment = all_data(all_data.ID == deployment_ids(i), :);
    
    % Identify indices where class is 'dest'
    dest_indices = find(strcmp(current_deployment.class, 'dest'));
    
    % Initialize arrays to store results
    start_idx = [];
    end_idx = [];
    mean_depths = [];
    std_depths = [];

    % Identify continuous sequences of 'dest'
    if ~isempty(dest_indices)
        % Find breaks in the sequence (where gap > 1)
        breaks = [1; diff(dest_indices) > 1]; % 1 where a new sequence starts
        start_idx = dest_indices(breaks == 1); % Start indices of sequences
        end_idx = [start_idx(2:end) - 1; dest_indices(end)]; % End indices
        
        % Compute mean and standard deviation of depth for each continuous 'dest' sequence
        mean_depths = zeros(size(start_idx)); % Ensure correct size
        std_depths = zeros(size(start_idx));  % Array for standard deviation
        for j = 1:length(start_idx)
            mean_depths(j) = mean(current_deployment.Depth(start_idx(j):end_idx(j)));
            std_depths(j) = std(current_deployment.Depth(start_idx(j):end_idx(j)));
        end
    end
    
    % Separate by Daylight classification ('Day' or 'Night')
    daylight_classification = current_deployment.Daylight; % 'Day' or 'Night'
    
    % Initialize containers for filtered results
    day_mean_depths = [];
    night_mean_depths = [];
    day_std_depths = [];
    night_std_depths = [];
    
    for k = 1:length(start_idx)
        % Check the corresponding 'Daylight' value for the current period
        daylights = daylight_classification(start_idx(k):end_idx(k));
        
        % For 'Day' classification
        if any(daylights == "Day")
            day_mean_depths = [day_mean_depths; mean_depths(k)];
            day_std_depths = [day_std_depths; std_depths(k)];
        end
        
        % For 'Night' classification (if needed)
        if any(daylights == "Night")
            night_mean_depths = [night_mean_depths; mean_depths(k)];
            night_std_depths = [night_std_depths; std_depths(k)];
        end
    end
    
    % Append the results of the current deployment to the arrays
    day_mean_depths_all = [day_mean_depths_all; day_mean_depths];
    night_mean_depths_all = [night_mean_depths_all; night_mean_depths];
    day_std_depths_all = [day_std_depths_all; day_std_depths];
    night_std_depths_all = [night_std_depths_all; night_std_depths];
end

% Calculate the overall mean and standard deviation across all deployments
overall_mean_day_depth = mean(day_mean_depths_all);
overall_std_day_depth = std(day_mean_depths_all);
overall_mean_night_depth = mean(night_mean_depths_all);
overall_std_night_depth = std(night_mean_depths_all);

% Display the overall results
disp(['Overall Mean Day Depth: ', num2str(overall_mean_day_depth)]);
disp(['Overall Standard Deviation of Day Depth: ', num2str(overall_std_day_depth)]);
disp(['Overall Mean Night Depth: ', num2str(overall_mean_night_depth)]);
disp(['Overall Standard Deviation of Night Depth: ', num2str(overall_std_night_depth)]);

%%
%Finding deep night feeeding dive kinematics 
%I guess I normalised the fluke rate so I'm just rerunning it here


% Replace NaN values in binary_peaks with 0 before performing the rolling sum
binary_peaks_no_nan = all_data.binary_peaks;
% Change every negative value to 0 and every positive value to 1
binary_peaks_no_nan(binary_peaks_no_nan < 0) = 0;  % Set negative values to 0
binary_peaks_no_nan(binary_peaks_no_nan > 0) = 1;   % Set positive values to 1


% Compute rolling sum of positive values in binary_peaks
all_data.fluking_event = movsum(binary_peaks_no_nan, [299 300]);



nightdeepfeed = all_data(all_data.MouthOpen == 1, :);
nightdeepfeed = nightdeepfeed(nightdeepfeed.Daylight == 'Night', :);
nightdeepfeed = nightdeepfeed(nightdeepfeed.Depth > 25, :);
disp('Speed');
disp(nanmean(nightdeepfeed.nnspeed));
disp('sd');
disp(nanstd(nightdeepfeed.nnspeed));

disp('Fluke Rate');
disp(nanmean(nightdeepfeed.fluking_event));
disp('sd');
disp(std(nightdeepfeed.fluking_event));
%%
night = all_data(strcmp(all_data.class, 'sur') & strcmp(all_data.Daylight, 'Night'), :);
sum(night.MouthOpen)/height(night);
%%
% Initialize arrays to store proportions for each ID and Hour
hourly_proportions_by_id = cell(24, 1);

% Loop over each hour (0 to 23)
for hour = 0:23
    % Initialize an array to hold the proportions for this hour
    proportions_for_hour = [];
    
    % Filter data for this hour
    hour_data = all_data(all_data.Hour == hour, :);
    
    % Loop through each unique ID
    unique_ids = unique(hour_data.ID);
    for id = unique_ids'
        % Filter rows for the current ID and hour
        id_data = hour_data(hour_data.ID == id, :);
        
        % Calculate the proportion of MouthOpen == 1 for this ID and hour
        mouth_open_data = id_data.MouthOpen == 1;
        proportion = sum(mouth_open_data) / height(id_data);
        
        % Append the proportion for this ID to the list for the current hour
        proportions_for_hour = [proportions_for_hour; proportion];
    end
    
    % Store the proportions for the current hour
    hourly_proportions_by_id{hour + 1} = proportions_for_hour;
end

% Calculate and display standard deviation for each hour
for hour = 0:23
    proportions_for_hour = hourly_proportions_by_id{hour + 1};
    
    % If there are multiple proportions, calculate the standard deviation
    if length(proportions_for_hour) > 1
        std_proportion = std(proportions_for_hour);
    else
        std_proportion = 0; % If only one value, the standard deviation is 0
    end
    
    disp(['Hour ', num2str(hour), ': Standard Deviation of Proportion of MouthOpen == 1 = ', num2str(std_proportion)]);
end

