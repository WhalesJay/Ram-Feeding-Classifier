
%%
%The purpose of this script is to apply a trained model, either the one
%provided in the package, one created using the dive_classify.mat script,
%or one created manually using the classifier toolbox in Matlab. 

%Section 1 - Data Loading
%This script requires the whale_names and *processed files from
%Dive_analyse.mat to be loaded. You can skip this section if you are
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
%Section 2 - Datenumber Conversion
%This section converts the datenumber into a more easily plottable format,
%the datenumber objects in matlab are digits which refer to a fixed
%reference data January 0, 0000 and are easier to keep stored like this
%until plotting. 

for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    whale_data = evalin('base', whale_name);
    whale_data.Datenum = datetime(whale_data.Datenum, 'ConvertFrom', 'datenum');
    whale_data.MouthOpen = NaN(height(whale_data), 1); % Adds an empty column (NaN for all rows) while we're here
    assignin('base', whale_name, whale_data);
end

clear("i", "whale_name")
disp("Section 2, Datenumber Conversion Complete");
%%
%Section 3 - Normalization
%Similar to other scripts in this package we will normalise certain metrics
%to enable inter-tag comparison 

% User-defined list of column names to search for
columns_to_search = {'Depth', 'Daylight', 'Hour', 'Datenum', 'nnspeed'};  % Modify this list as needed

for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name (assuming whale_names is a cell array)
    
    % Access the whale_data object dynamically from the base workspace
    whale_data = evalin('base', whale_name); 
    
    % Save 'speed' for filtration estimates later
    whale_data.nnspeed = whale_data.speed; 
    
    % Find columns to normalize by checking if they exist in whale_data
    columns_to_normalize = ismember(whale_data.Properties.VariableNames, columns_to_search);
    
    % Loop through columns and normalize if they are in the list
    for j = 1:width(whale_data)
        if ~columns_to_normalize(j)
            % Normalize the column by subtracting mean and dividing by standard deviation
            whale_data{:, j} = (whale_data{:, j} - mean(whale_data{:, j}, 'omitnan')) / std(whale_data{:, j}, 'omitnan');
        end
    end
    
    % Assign the modified whale_data back to the base workspace
    assignin('base', whale_name, whale_data);
end

clear("i", "j", "columns_to_normalize", "whale_data", "whale_name", "columns_to_search");
disp("Section 3, Normalisation - Complete");
%%
%Section 4 - Application of Trained model
%Run and rerun this section if you are trialling multiple models being sure
%to either replace the model name or replace the model item in the
%workspace
% Extract data from whale_data
% Loop through each whale name in whale_names

load("trainedModel.mat"); %Replace this with your own or another trained model
for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    whale_data = evalin('base', whale_name);
    %whale_data = renamevars(whale_data, 'hour', 'Hour');
    [mouthOpen, scores] = medSVM.predictFcn(whale_data); %change this model name if exporting direct to workspace from classification toolbox
    whale_data.MouthOpen = mouthOpen;

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
        % Create the file name
    file_name = [whale_name, '_SVM_Classified.mat'];

    % Save the file
    save(file_name, 'whale_data');
end    
nn_data = table(); % Initialize an empty table to hold data from all whale names

% Loop through each whale_name object and concatenate the specified columns
for i = 1:length(whale_names)
    whale_name = whale_names(i).name; % Get the current whale name
    % Access the whale_data object dynamically
    current_data = evalin('base', whale_name);
    current_data.WhaleName = repmat(string(whale_name), height(current_data), 1);
    nn_data = [nn_data; current_data]; % Append to all_data table
end

clear("whale_data", "i", "Date", "Depth", "MouthOpen");
disp("Section 4, Model Fitting & Plotting - Complete");
%%
%Section 5 - Model Evaluation
%This section uses the training data to calculate various evaluation
%metrics and can be skipped if you do not need to do any further model
%evaluation 
models = {bagged, bilayeredNN, fineSVM, mediumNN, medSVM, wideNN, trilayeredNN};
modelNames = {'bagged', 'bilayeredNN', 'fineSVM', 'mediumNN', 'medSVM', 'wideNN', 'trilayeredNN'};
metrics = table('Size', [length(models), 5], ...
                'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'Model', 'Accuracy', 'Precision', 'Recall', 'F1Score'});
load('normData.mat')


for i = 1:length(models)
    model = models{i};
    modelName = modelNames{i};
    
    % Determine how to extract X and Y based on model type
    if contains(modelName, 'NN')  % Neural Networks
        XTrain = model.ClassificationNeuralNetwork.X;
        YTrain = model.ClassificationNeuralNetwork.Y;
    elseif contains(modelName, 'SVM')  % SVM
        XTrain = model.ClassificationSVM.X;
        YTrain = model.ClassificationSVM.Y;
    else  % Default to bagged (or others not matching above)
        XTrain = model.ClassificationEnsemble.X;
        YTrain = model.ClassificationEnsemble.Y;
    end
       % Add blank 'Date' column if it does not exist
    if ~ismember('Date', XTrain.Properties.VariableNames)
        XTrain.Date = NaN(height(XTrain), 1);  % Add blank 'Date' column with NaN values
    end
    % Get predictions using the model's predict function
    YPred = model.predictFcn(XTrain);
    
    % Compute confusion matrix and metrics
    cm = confusionmat(YTrain, YPred);
    TP = diag(cm);  % True Positives
    FP = sum(cm, 1)' - TP;  % False Positives
    FN = sum(cm, 2) - TP;  % False Negatives
    
    % Handle cases where precision or recall might be NaN (e.g., division by zero)
    precision = mean(TP ./ (TP + FP + eps));  % Add eps to avoid NaN
    recall = mean(TP ./ (TP + FN + eps));
    f1score = mean(2 * (precision .* recall) ./ (precision + recall + eps));
    accuracy = sum(TP) / sum(cm(:));
    
    %metrics{i, 1} = modelName;  % Model name (string)
    metrics{i, 2} = accuracy;   % Accuracy (double)
    metrics{i, 3} = precision;  % Precision (double)
    metrics{i, 4} = recall;     % Recall (double)
    metrics{i, 5} = f1score;    % F1 Score (double)
end

disp(metrics);

disp("Section 5 - Model Evaluation Complete")
%%
%Section 6 - Mouth Open Characteristics 

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


all_data.MouthOpen = categorical(all_data.MouthOpen, [0 1], {'Closed', 'Open'});

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
            v = violinplot(all_data.(column_name), all_data.MouthOpen);
            
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
disp("Section 6, Boxplots - Complete")



%%
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

        % Compute histogram counts and bin edges for depth, normalize by percentage
        [counts, binEdges] = histcounts(currentDepth, edges);
        counts = (counts / totalRows) * 100;  % Convert frequency to percentage
        bar(binEdges(1:end-1), counts, 'FaceColor', '#000100', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        hold on;
        
        % Compute histogram counts and bin edges for speed, normalize by percentage
        [counts, binEdges] = histcounts(currentSpeed, edges);
        counts = (counts / totalRows) * 100;  % Convert frequency to percentage
        bar(binEdges(1:end-1), counts, 'FaceColor', '#FB3640', 'EdgeColor', 'none', 'FaceAlpha', 0.8);

        % Compute histogram counts and bin edges for fluking event, normalize by percentage
        [counts, binEdges] = histcounts(currentFluking, edges);
        counts = (counts / totalRows) * 100;  % Convert frequency to percentage
        bar(binEdges(1:end-1), counts, 'FaceColor', '#48ACF0', 'EdgeColor', 'none', 'FaceAlpha', 0.6);


        % Add labels
        xlabel('Normalized Value', 'FontSize', 10);
        ylabel('Percentage (%)', 'FontSize', 10);
        
        % Adjust x-axis limits for uniformity across all subplots
        xlim([0, 1]);
        ylim([0,18]);
        
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
text(-0.3, 0.5, 'Day', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
subplot(3, 3, 4);
text(-0.3, 0.5, 'Night', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
subplot(3, 3, 7);
text(-0.3, 0.5, 'Twilight', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% Add a single overall legend outside of the subplots (to the right)
legend({'Speed', 'Fluke Rate', 'Depth'}, 'Location', 'BestOutside', 'FontSize', 10, 'Box', 'off', 'Orientation', 'vertical');

% Adjust the overall layout and add a super title
set(gcf, 'Color', 'w');  % Set figure background to white for clarity

% Manually adjust legend position to be outside of the plot
legendPosition = [0.92, 0.5, 0.05, 0.1];  % Adjust this to move the legend to the desired location
legend('Position', legendPosition);
%% Dive Patterns
%This script visualises dive depth through a 24h cycle
%This requires the all_data output from the Dive_visualise.m script 

hour_bins = 0:23;  % 24 hours of the day
depth_bins = linspace(0, 80, 20);  % 30 depth bins between 0 and 80 meters

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
% Adjust layout for clear separation and visibility
%%
%This is just a work around because I am struggling to reverse the
%normalisation process, so I saved the data before normalising, now I am
%just adding back in the mouth open prediction and then calculating the
%mean and std for each variable 
% Add MouthOpen column to nn_data and convert to binary (1 for 'Open', 0 for 'Closed')
nn_data.MouthOpen = all_data.MouthOpen;
nn_data.MouthOpen = double(nn_data.MouthOpen == "Open");

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

% Display results
disp('Average Results by Daylight and MouthOpen:');
disp(averageResults);

disp('Standard Deviation Results by Daylight and MouthOpen:');
disp(stdResults);



%%
%Making some summary stats

%First is the percentage of time with mouth open when in the top 5m of the
%water column 

% Validate that the required columns exist in nn_data
requiredColumns = {'MouthOpen', 'Depth', 'Daylight'};
if all(ismember(requiredColumns, nn_data.Properties.VariableNames))
    % Filter rows where Depth < 5
    depthFilter = nn_data.Depth > 20;

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
