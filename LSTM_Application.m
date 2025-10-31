%This script allows the user to apply the results of either their own or
%the included pre-trained LSTM model to unannotated data. 


%Section 1 - Data Initiation 
%This section loads up all the deployments required for classification
%using the previously created whale_names.mat and processed files from
%Dive_analyse.mat. You can skip this section if you are
%running immediately after running the Dive_analyse.mat script without
%clearing the workspace

% Part 1: Load and Store whale_data Structures
whale_names = dir('*processed.mat');

% Loop through the list and filter out files with underscores in their names
for i = 1:length(whale_names)
    % Check if the file name does not contain an underscore
    if isempty(strfind(whale_names(i).name, '_'))
        % Load the file if it meets the criteria
        data = load(whale_names(i).name);
        
        % Create a variable name without the .mat extension
        var_name = erase(whale_names(i).name, '.mat');
        
        % Assign the whale_data structure to the base workspace
        assignin('base', var_name, data.whale_data);
    end
end

clear("whale_names", "var_name", "i", "data");

%Load in the whale names list
load('whale_names.mat'); 

disp('Section 1, File Loading - Complete');
%%
%%EDIT REQ
%Section 2 - Variable Selection
%This section removes the variables which were not used in the LSTM,
%leaving only the accelerometry, gyroscope, roll, heading, speed, jerk and
%derived pitch variables (body orientation and fluking signal). 

%The loop then creates sequences for prediction and then uses the trained
%model to remove all 

sequenceLength = 300;  % Length of each sequence
overlapStep = 1;       % Overlap by 1 row
load('mouthopenLSTM.mat');  % Load the trained LSTM model
numRuns = (sequenceLength/overlapStep); %Number of times the prediction will rerun
startPoints = (0:numRuns-1) * (overlapStep) + 1;

for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    whale_data = evalin('base', whale_name);
    % Keep 'Datenum' and 'Depth' columns in whale_data but exclude them from classification
    %valid_rows = all(isfinite(table2array(whale_data)), 2); % Logical index for valid rows
    %whale_data = whale_data(valid_rows, :); % Keep only valid rows
    %Remove rows with NAs (Usually speed or jerk at start of deployment)
    numericCols = varfun(@isnumeric, whale_data, 'OutputFormat', 'uniform');
    numericData = table2array(whale_data(:, numericCols));
    valid_rows = all(isfinite(numericData), 2);
    whale_data = whale_data(valid_rows, :);
    numRows = height(whale_data);
    predictionsMatrix = NaN(numRows, numRuns); % Rows: data, Columns: start points
    probabilitiesMatrix = NaN(numRows, numRuns);
        for g = 1:numRuns
% Segment the data into sequences of 300 rows each
        numSequences = floor((height(whale_data)-startPoints(g)) / sequenceLength); % Number of complete sequences
        inputData = cell(numSequences, 1);  % Cell array for storing input sequences
        labels = cell(numSequences, 1);    % Cell array for storing the corresponding labels
        predicted_labels = NaN(height(whale_data), 1);  % Initialize column for predictions (NaN for all rows)
            for j = 1:numSequences
            % Get the start and end indices for the current sequence
            startIdx = ((((j-1)*sequenceLength)+startPoints(g)));
            endIdx = startIdx + sequenceLength-1;
            % Extract the 18 features (excluding 'Datenum' and 'Depth') for the sequence
            % Each sequence is a 300x14 matrix, but we want it as 14x300 (transposed)
            inputData{j} = table2array(whale_data(startIdx:endIdx, [3, 4, 5, 16, 15, 6, 7, 8, 9, 12, 11, 13, 18,19]))';
            end
        % Apply the trained LSTM model to each sequence
        [predictions, probabilities] = classify(net, inputData);
        maxPredictRow = startPoints(g) + numSequences * 300;
        %predictionsMatrix([startPoints(g):maxPredictRow],g) = predictions;
        %probabilitesMatrix([startPoints(g):maxPredictRow],g) = probabilities(:,1);
        % Get the number of rows in whale_data
        numRows = height(whale_data);
        % Start index for assigning values
        v = startPoints(g);
% Loop through the predictions and assign in chunks of 300
            for l = 1:height(predictions)
            b = v + 299;  % End index for current block
            block = str2double(string(predictions(l)));
            predictionsMatrix(v:b, (g)) = block;  % Assign the prediction to the chunk
            block = str2double(string(probabilities(l)));
            probabilitiesMatrix(v:b, (g)) = block;  % Assign the prediction to the chunk
            v = b + 1;  % Update start index for next block
            end
        end
% Set the MouthOpen column to NaNs
    whale_data.MouthOpen = mode(predictionsMatrix, 2);
    whale_data.MouthOpenProb = nanmean(probabilitiesMatrix, 2);
    
    % Store the predictions in the base workspace
    assignin('base', [whale_name, '_LSTM_Classified'], whale_data);

    % Create the file name
    file_name = [whale_name, '_LSTM_Classified.mat'];

    % Save the file
    save(file_name, 'inputData', 'labels', 'predictions', 'whale_data');
    current_time = datetime('now', 'Format', 'HH:mm:ss');
    
    % Display the whale_name along with the timestamp
    disp([whale_name, ' classified at ', char(current_time)]);
end

whale_names = whos("*LSTM_Classified");
% Extract the names of these variables into a cell array
variables_to_keep = {whale_names.name};

% Add 'whale_names' itself to the list of variables to keep
variables_to_keep = [variables_to_keep, 'whale_names'];

% Clear all other variables except those in 'variables_to_keep'
clearvars('-except', variables_to_keep{:});

disp("Section 2 -Predictions Made")
%%
%Section 3 - Plotting Predictions
%This section is for plotting the outputs of the classification
%predicitions from the above section. This section can be changed to plot
%probabilities on line 138


for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    whale_data = evalin('base', whale_name);

    mouthOpen = whale_data.MouthOpen; %Change to whale_data.MouthOpenProb to plot probabilities

    %whale_data.Datenum = datetime(whale_data.Datenum, 'ConvertFrom', 'datenum');
    % Extract data from whale_data
    Date = whale_data.Datenum;  % Assuming Date is a datetime or numeric vector
    Depth = whale_data.Depth;  % Depth values

    % Create a figure for each whale
    figure_handle = figure;

    % Create a scatter plot with Date on x-axis and Depth on y-axis
    scatter(Date, Depth, 36, mouthOpen, 'filled');  % 36 is the marker size

    % Invert the y-axis (to represent increasing depth downwards)
    set(gca, 'YDir', 'reverse');

    % Add colorbar to indicate the mapping of mouthOpen values
    colorbar;
    colormap(jet);  % You can choose other colormaps like 'hot', 'parula', etc.

    % Add labels and title
    xlabel('Date');
    ylabel('Depth');
    title(sprintf('Whale Depth vs. Date for %s with Mouth Open Indicator', whale_name));
    saveas(figure_handle, strcat(whale_name, '_classified.png')); % Save as PNG format
    assignin('base', whale_name, whale_data);
end    

clear("whale_data", "i", "Date", "Depth", "MouthOpen");
disp("Section 3, Model Plotting - Complete");
%% 
%Section 4 - Mouth Open Characteristics 

%This section first creates a large data object similar to the one used in
%the Dive_pattern.mat script, then creates boxplots comparing mouth open
%and closed for each of the predictors based on the the predicted value of
%MouthOpen

all_data = table(); % Initialize an empty table to hold data from all whale names

% Loop through each whale_name object and concatenate the specified columns
for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    current_data = evalin('base', whale_name);
    current_data.WhaleName = repmat(string(whale_name), height(current_data), 1);
    all_data = [all_data; current_data]; % Append to all_data table
end

clear("current_data", "whale_name", "i");


%all_data.MouthOpen = categorical(all_data.MouthOpen, [0 1], {'Closed', 'Open'});

%Index the columns to plot
exclude_columns = find(ismember(all_data.Properties.VariableNames, {'Date', ['Datenum",' ...
    'binary_peaks', 'Daylight', 'WhaleName', 'Hour']}));
% Define a softer color palette (e.g., pastel colors)
colors = [0.7 0.8 0.9; 0.9 0.7 0.8]; % Example of soft pastel colors


% Loop through each column in the table, excluding specified columns
for i = 1:width(all_data)
    if ~ismember(i, exclude_columns)
        column_name = all_data.Properties.VariableNames{i};
        
        % Skip the 'MouthOpen' column itself
        if ~strcmp(column_name, 'MouthOpen') && isnumeric(all_data.(column_name))
            % Create a figure for each column
            figure;
            
            % Create the violin plot with the specified colors
            valid_idx = ~isnan(all_data.(column_name)) & ~isnan(all_data.MouthOpen);
                if any(valid_idx)
                    v = violinplot(all_data.(column_name)(valid_idx), all_data.MouthOpen(valid_idx));
                else
                    error('No valid data to plot.');
                end
            
            % Title and labeling
            title(['Violin plot of ', column_name, ' by MouthOpen']);
            xlabel('MouthOpen');
            ylabel(column_name);
            
            % Save the plot with the title as the filename
            saveas(gcf, ['Violin_plot_', column_name, '.png']);
        end
    end
end

clear("column_name", "exclude_columns", "figure_handle", "i", "j", "colors", "v");
disp("Section 4, Boxplots - Complete")

%%
%This section creates our key variable plots
% Ensure the variables are in the correct format
speed = all_data.speed;
fluking_event = all_data.fluking_event;
depth = all_data.Depth;

% Normalize each variable by its max value to ensure shared scale
speed_norm = speed / max(speed);
fluking_event_norm = fluking_event / max(fluking_event);
depth_norm = depth / max(depth);

% Ensure that 'Daylight' and 'MouthOpen' are categorical
daylight = categorical(all_data.Daylight);  % Daylight with 3 categories: 'Day', 'Twilight', 'Night'
mouthOpen = categorical(all_data.MouthOpen); % MouthOpen with 2 categories: 0 or 1

% Define consistent color palette for each variable
speedColor = [200/255, 30/255, 34/255];  % Green
flukingColor = [30/255, 144/255, 255/255];  % Blue
depthColor = [169/255, 169/255, 169/255];  % Grey

% Get the categories (they are now cell arrays)
daylightCategories = categories(daylight);
mouthOpenCategories = categories(mouthOpen);

% Set up a figure with subplots for each combination
figure('Position', [100, 100, 1200, 900]);

% Loop through each combination of Daylight and MouthOpen
for i = 1:length(daylightCategories)
    for j = 1:length(mouthOpenCategories)
        
        % Select data for the current combination of Daylight and MouthOpen
        currentDaylightIdx = daylight == daylightCategories{i};
        currentMouthOpenIdx = mouthOpen == mouthOpenCategories{j};
        currentIdx = currentDaylightIdx & currentMouthOpenIdx;  % Combine the two conditions
        
        % Get the corresponding data for each variable
        currentSpeed = speed_norm(currentIdx);
        currentFluking = fluking_event_norm(currentIdx);
        currentDepth = depth_norm(currentIdx);
        
        % Create a subplot for this combination (3 rows and 3 columns of histograms)
        subplot(3, 3, (i-1)*3 + j);
        
        % Define the bins and the edges (smaller step size for more bins)
        edges = 0:0.0125:1;  % 0 to 1 with 0.0125 steps (4x smaller than 0.05)

        % Get the total number of rows in the dataset (to calculate percentage)
        totalRows = sum(currentIdx);

        % Compute histogram counts and bin edges for speed, normalize by percentage
        [counts, binEdges] = histcounts(currentSpeed, edges);
        counts = (counts / totalRows) * 100;  % Convert frequency to percentage
        bar(binEdges(1:end-1), counts, 'FaceColor', speedColor, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        hold on;

        % Compute histogram counts and bin edges for fluking event, normalize by percentage
        [counts, binEdges] = histcounts(currentFluking, edges);
        counts = (counts / totalRows) * 100;  % Convert frequency to percentage
        bar(binEdges(1:end-1), counts, 'FaceColor', flukingColor, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

        % Compute histogram counts and bin edges for depth, normalize by percentage
        [counts, binEdges] = histcounts(currentDepth, edges);
        counts = (counts / totalRows) * 100;  % Convert frequency to percentage
        bar(binEdges(1:end-1), counts, 'FaceColor', depthColor, 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Add labels
        xlabel('Normalized Value', 'FontSize', 10);
        ylabel('Percentage (%)', 'FontSize', 10);
        
        % Adjust x-axis limits for uniformity across all subplots
        xlim([0, 1]);
        ylim([0,25]);
        
        % Turn off scientific notation on the y-axis and set the format for the percentages
        ax = gca;
        ax.YAxis.Exponent = 0;  % Disable scientific notation
        ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), ax.YTick, 'UniformOutput', false);  % Format as percentage
        
        hold off;
    end
end
% Add titles to the columns (Mouth Closed / Mouth Open)
subplot(3, 3, 1);
title('Mouth Closed', 'FontSize', 14, 'FontWeight', 'bold');
subplot(3, 3, 2);
title('Mouth Open', 'FontSize', 14, 'FontWeight', 'bold');

% Add titles to the rows (Day, Night, Twilight)
subplot(3, 3, 1);
text(-0.3, 0.5, 'Day', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
subplot(3, 3, 4);
text(-0.3, 0.5, 'Night', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
subplot(3, 3, 7);
text(-0.3, 0.5, 'Twilight', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% Add a single overall legend outside of the subplots (to the right)
legend({'Speed', 'Fluking Event', 'Depth'}, 'Location', 'BestOutside', 'FontSize', 10, 'Box', 'off', 'Orientation', 'vertical');

% Adjust the overall layout and add a super title
sgtitle('Frequency Distribution of Speed, Fluking Event, and Depth by Daylight and Mouth Open Classification', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');  % Set figure background to white for clarity

% Manually adjust legend position to be outside of the plot
legendPosition = [0.92, 0.5, 0.05, 0.1];  % Adjust this to move the legend to the desired location
legend('Position', legendPosition);

%%
%% Dive Patterns
%This script visualises dive depth through a 24h cycle
%This requires the all_data output from the Dive_visualise.m script 

hour_bins = 0:23;  % 24 hours of the day
depth_bins = linspace(0, 80, 30);  % 30 depth bins between 0 and 80 meters

% Ensure MouthOpen is categorical for correct processing
mouthOpen = categorical(all_data.MouthOpen);  % Assuming 'Open' and 'Closed' are categories

% Filter the data based on MouthOpen categories
mouthClosedData = all_data(mouthOpen == '0', :);
mouthOpenData = all_data(mouthOpen == '1', :);

% Step 5: Bin the data and calculate density for Mouth Closed
[counts_closed, hour_edges_closed, depth_edges_closed] = histcounts2(mouthClosedData.Hour, mouthClosedData.Depth, hour_bins, depth_bins);

% Step 6: Normalize counts for Mouth Closed by points per hour bin
points_per_hour_closed = sum(counts_closed, 2);  % Sum across depth bins for each hour
points_per_hour_closed(points_per_hour_closed == 0) = NaN;  % Avoid division by zero
normalized_counts_closed = counts_closed ./ points_per_hour_closed;  % Normalize

% Step 7: Bin the data and calculate density for Mouth Open
[counts_open, hour_edges_open, depth_edges_open] = histcounts2(mouthOpenData.Hour, mouthOpenData.Depth, hour_bins, depth_bins);

% Step 8: Normalize counts for Mouth Open by points per hour bin
points_per_hour_open = sum(counts_open, 2);  % Sum across depth bins for each hour
points_per_hour_open(points_per_hour_open == 0) = NaN;  % Avoid division by zero
normalized_counts_open = counts_open ./ points_per_hour_open;  % Normalize

% Optional: Replace NaNs with zeros for visualization purposes
normalized_counts_closed(isnan(normalized_counts_closed)) = 0;
normalized_counts_open(isnan(normalized_counts_open)) = 0;

% Now, normalized_counts_closed and normalized_counts_open contain density normalized to points per hour bin for each mouth state

% Step 7: Create the figure with two subplots for the heatmaps
figure;

% Subplot 1: Heatmap for Mouth Closed
subplot(1, 2, 1);  % 1 row, 2 columns, first subplot
imagesc(hour_edges_closed(1:end-1), depth_edges_closed(1:end-1), log(counts_closed)');  % Transpose counts to align with axes
colorbar off;  % Turn off the colorbar for the first subplot
xlabel('Hour of Day','FontSize', 20);
ylabel('Depth (m)','FontSize', 20);
set(gca, 'FontSize', 14);
title('Mouth Closed','FontSize', 30);
colormap('parula');  % Use a color map for better visualization
set(gca, 'YDir', 'reverse');  % Reverse the y-axis to have depth increasing downward


% Subplot 2: Heatmap for Mouth Open
subplot(1, 2, 2);  % 1 row, 2 columns, second subplot
imagesc(hour_edges_open(1:end-1), depth_edges_open(1:end-1), log(counts_open)');  % Transpose counts to align with axes
cb = colorbar;  % Show colorbar for density
xlabel('Hour of Day','FontSize', 20);
set(gca, 'FontSize', 14);
set(gca, 'YTickLabel', []);
title('Mouth Open', 'FontSize', 30);
colormap('parula');  % Use a color map for better visualization
set(gca, 'YDir', 'reverse');  % Reverse the y-axis to have depth increasing downward

cb.Title.String = 'Log Density of time spent in 2m depth bin per hour';% This sets the title of the colorbar.
cb.Title.FontSize = 15; % Optionally sets the font size of the colorbar title.
cb.Title.Rotation = 270;  % Rotate title by 90 degrees
cb.Title.Position = [30, 175, 0];  % Adjust title position to be next to the colorbar
%colorbar('off');


%%
%%
%Making summary stats for each variable for each daylight and mouth
%classification
nn_data = all_data;
% Retain only numeric columns, excluding 'MouthOpen'
numericData = nn_data(:, varfun(@isnumeric, nn_data, 'OutputFormat', 'uniform'));
numericData = removevars(numericData, 'MouthOpen');  % Retain numeric columns but remove MouthOpen

% Ensure Daylight is handled as a categorical variable (if it's a string or cell)
if iscell(nn_data.Daylight) || ischar(nn_data.Daylight)
    nn_data.Daylight = categorical(nn_data.Daylight);
end

% Extract unique categories in Daylight
uniqueDaylight = unique(nn_data.Daylight);

% Initialize tables to store results (dynamic allocation)
averageResults = table();
stdResults = table();

% Loop through each Daylight category
for i = 1:length(uniqueDaylight)
    % Filter data for the current Daylight category
    currentDaylight = uniqueDaylight(i);
    daylightFilter = nn_data.Daylight == currentDaylight;

    % Further split data based on MouthOpen values
    dataDaylight = numericData(daylightFilter, :);
    mouthOpen1 = nn_data.MouthOpen(daylightFilter) == 1;
    mouthOpen0 = nn_data.MouthOpen(daylightFilter) == 0;

    dataMouthOpen1 = dataDaylight{mouthOpen1, :};
    dataMouthOpen0 = dataDaylight{mouthOpen0, :};

    % Calculate mean and standard deviation for the current Daylight category
    meanMouthOpen1 = mean(dataMouthOpen1, 'omitnan');
    meanMouthOpen0 = mean(dataMouthOpen0, 'omitnan');
    stdMouthOpen1 = std(dataMouthOpen1, 'omitnan');
    stdMouthOpen0 = std(dataMouthOpen0, 'omitnan');

    % Append results to the tables for each row
    averageResults = [averageResults; table(repmat(string(currentDaylight), size(meanMouthOpen1, 2), 1), meanMouthOpen1', meanMouthOpen0', ...
        'VariableNames', {'Daylight', 'MouthOpen1', 'MouthOpen0'})];

    stdResults = [stdResults; table(repmat(string(currentDaylight), size(stdMouthOpen1, 2), 1), stdMouthOpen1', stdMouthOpen0', ...
        'VariableNames', {'Daylight', 'MouthOpen1', 'MouthOpen0'})];
end



%%
%Making some summary stats

%First is the percentage of time with mouth open when in the top 5m of the
%water column 

% Validate that the required columns exist in nn_data
requiredColumns = {'MouthOpen', 'Depth', 'Daylight'};
if all(ismember(requiredColumns, nn_data.Properties.VariableNames))
    % Filter rows where Depth < 5
    depthFilter = nn_data.Depth <5;

    % Group by Daylight categories
    uniqueDaylight = unique(nn_data.Daylight); % Get unique Daylight categories

    % Initialize results table
    daylightCategories = cellstr(uniqueDaylight); % Convert categories to cell array of strings
    percentages = zeros(size(uniqueDaylight)); % Preallocate array for percentages

    % Loop through each Daylight category
    for i = 1:length(uniqueDaylight)
        % Filter rows for the current Daylight category
        daylightFilter = nn_data.Daylight == uniqueDaylight(i);

        % Rows satisfying both Depth < 5 and the current Daylight category
        combinedFilter = depthFilter & daylightFilter;

        % Total rows in the current Daylight category with Depth < 5
        totalRows = sum(daylightFilter & depthFilter);

        % Rows with MouthOpen = 1 in the current category and Depth < 5
        mouthOpenRows = sum(nn_data.MouthOpen == 1 & combinedFilter);

        % Calculate percentage
        if totalRows > 0
            percentages(i) = (mouthOpenRows / totalRows) * 100;
        else
            percentages(i) = NaN; % Handle cases with no data
        end
    end

    % Create results table
    results = table(daylightCategories, percentages, ...
        'VariableNames', {'Daylight', 'Percentage_MouthOpen_1'});

    % Display results
    disp('Percentage of rows with MouthOpen = 1 and Depth < 5 by Daylight:');
    disp(results);
else
    error('One or more required columns (MouthOpen, Depth, Daylight) are missing in nn_data.');
end
%%
%Finding average feeding duration and depth of feeding
% Define daylight categories manually
daylight_labels = ["avg_depth_day", "avg_depth_night", "avg_depth_twilight"];
prop_labels = ["prop_mouthOpen_day", "prop_mouthOpen_night", "prop_mouthOpen_twilight"];

% Define column names
col_names = ["whale_name", "avg_depth_overall", daylight_labels, "prop_mouthOpen_overall", prop_labels];

% Initialize an empty results table with correct column types
results_table = table('Size', [0, numel(col_names)], ...
                      'VariableTypes', ["string", repmat("double", 1, numel(col_names) - 1)], ...
                      'VariableNames', col_names);

for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    whale_data = evalin('base', whale_name); % Access whale_data dynamically

    % Convert whale_name to string for table storage
    whale_name_str = string(whale_name);
    
    % Overall average Depth when MouthOpen == 1
    avg_depth_overall = mean(whale_data.Depth(whale_data.MouthOpen == 1), 'omitnan');
    
    % Proportion of rows where MouthOpen == 1 (overall)
    prop_mouthOpen_overall = mean(whale_data.MouthOpen == 1, 'omitnan');
    
    % Compute values for each daylight category
    avg_depth_day = mean(whale_data.Depth(whale_data.MouthOpen == 1 & whale_data.Daylight == "Day"), 'omitnan');
    avg_depth_night = mean(whale_data.Depth(whale_data.MouthOpen == 1 & whale_data.Daylight == "Night"), 'omitnan');
    avg_depth_twilight = mean(whale_data.Depth(whale_data.MouthOpen == 1 & whale_data.Daylight == "Twilight"), 'omitnan');

    prop_mouthOpen_day = mean(whale_data.MouthOpen(whale_data.Daylight == "Day") == 1, 'omitnan');
    prop_mouthOpen_night = mean(whale_data.MouthOpen(whale_data.Daylight == "Night") == 1, 'omitnan');
    prop_mouthOpen_twilight = mean(whale_data.MouthOpen(whale_data.Daylight == "Twilight") == 1, 'omitnan');

    % Create new row for the table
    new_row = table(whale_name_str, avg_depth_overall, avg_depth_day, avg_depth_night, avg_depth_twilight, ...
                    prop_mouthOpen_overall, prop_mouthOpen_day, prop_mouthOpen_night, prop_mouthOpen_twilight, ...
                    'VariableNames', col_names);

    % Append to results table
    results_table = [results_table; new_row];
end

% Display the results table
disp(results_table);
