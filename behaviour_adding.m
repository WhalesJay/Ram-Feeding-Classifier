%Section 1 - Loading Behaviour Data
%Once the behavioural audit has been completed from BORIS, export the
%behaviours as aggregated events, selecting yes to add them all to one
%sheet. Before loading the tsv from BORIS open the file in Excel and alter the
%"Media File Name" column so that it reads only as the name of the video
% (e.g. split by "/" and keep only the filename). Then save as a csv
% in the working directory.

% Define the CSV file to read, change this!
filename = 'training.csv';

% Read the CSV file, treating unreadable cells as missing values (NaN)
opts = detectImportOptions(filename);

% Read the table with the modified options
behav = readtable(filename, opts);
clear("filename", "opts");

disp("Section 1, Behaviour Data Loaded");

%% Section 2 - Set fixed FPS and extract video start time (vidst) from MediaFileName
% Frame rate is fixed at 25 fps for all videos.
fixed_fps = 30;
behav.NewFPS = repmat(fixed_fps, height(behav), 1);

% Initialize vidst (video start datetime) as NaT
behav.vidst = NaT(height(behav), 1);

% Regex patterns we will check:
% 1) yyyyMMdd-HHmmss  -> '20250930-115812'
% 2) yyyyMMddHHmmss   -> '20250930115812'
pat1 = '(\d{8}-\d{6})';
pat2 = '(\d{14})';

for i = 1:height(behav)
    fname = behav.MediaFileName{i};
    if isempty(fname), continue; end
    
    % try pattern 1: YYYYMMDD-HHMMSS
    tkn = regexp(fname, pat1, 'match', 'once');
    if ~isempty(tkn)
        try
            behav.vidst(i) = datetime(tkn, 'InputFormat', 'yyyyMMdd-HHmmss');
            continue; % proceed next row
        catch
            % fall through to try other patterns
        end
    end
    
    % try pattern 2: YYYYMMDDHHMMSS (no hyphen)
    tkn = regexp(fname, pat2, 'match', 'once');
    if ~isempty(tkn)
        try
            behav.vidst(i) = datetime(tkn, 'InputFormat', 'yyyyMMddHHmmss');
            continue;
        catch
            % leave as NaT
        end
    end
    
    % If no timestamp found, leave behav.vidst(i) as NaT
end

% Optional: try to extract a sequence number from the filename if behav.seq
% does not already exist (useful for Section 6's "next video" logic).
% This attempts to read the last contiguous run of digits before the file
% extension or trailing dots (e.g. '...-000-00004..mp4' -> 4).
if ~ismember('seq', behav.Properties.VariableNames)
    behav.seq = NaN(height(behav),1);
    for i = 1:height(behav)
        fname = behav.MediaFileName{i};
        if isempty(fname), continue; end
        % Capture last run of digits before the extension or trailing dots
        seq_str = regexp(fname, '(\d+)(?=(\.[^.]*$)|\.\.)', 'match', 'once');
        if ~isempty(seq_str)
            % convert to number (strip leading zeros)
            seqnum = str2double(seq_str);
            if ~isnan(seqnum)
                behav.seq(i) = seqnum;
            end
        end
    end
end

disp("Section 2, Fixed FPS (30) Applied and vidst extracted from filenames");

%% Section 3 - Compute frame-based timing using fixed FPS
% Convert frame indices to seconds using the fixed fps and add to vidst.

behav.fstart = behav.ImageIndexStart ./ fixed_fps;
behav.fstop  = behav.ImageIndexStop  ./ fixed_fps;

% Initialize tstart and tstop as datetime arrays (already done earlier but
% re-init here to be safe)
behav.tstart = NaT(height(behav), 1);
behav.tstop  = NaT(height(behav), 1);

for i = 1:height(behav)
    if ~isnat(behav.vidst(i))
        behav.tstart(i) = behav.vidst(i) + seconds(behav.fstart(i));
        behav.tstop(i)  = behav.vidst(i) + seconds(behav.fstop(i));
    else
        % vidst missing -> tstart/tstop remain NaT (warn optionally)
        % fprintf('Note: vidst missing for %s (row %d)\n', behav.MediaFileName{i}, i);
    end
end

disp("Section 3, Frame Timing Computed from vidst and fixed FPS");
%% Section 4 - Load Fully Audited PRH Files
% The following sections load and finish the processing for the PRH files

processed_files = dir('*processed.mat');

for i = 1:length(processed_files)
    file_name = processed_files(i).name;
    mat_data = load(file_name);
    [~, whale_name, ~] = fileparts(file_name);

    if isfield(mat_data, 'whale_data')
        assignin('base', whale_name, mat_data.whale_data);
    else
        warning(['whale_data not found in file: ', file_name]);
    end
end

whale_names = evalin('base', 'whos');
whale_names = whale_names(startsWith({whale_names.name}, 'mn'));

clear("processed_files", "file_name", "i", "mat_data", "whale_name");
disp("Section 4, Audited PRH Files Loaded");

%% Section 5 - Convert Datetime to plottable numbers
for i = 1:length(whale_names)
    whale_name = whale_names(i).name;
    whale_data = eval(whale_name);
    whale_data.Datenum = datetime(whale_data.Datenum, 'ConvertFrom', 'datenum');
    assignin('base', whale_name, whale_data);
end

clear('i', "whale_data", "whale_name");
disp("Section 5 - Datetime Converted");

%% Section 6 - Extract unique behaviors
unique_behaviors = unique(behav.Behavior);

observation_mapping = containers.Map( ...
    {'mn42_30092025', 'rw44_240713', 'rw45_07072024', 'rw69_11072024'}, ...
    {'mn3009202542processed', 'rw24071344processed', 'rw24070745processed', 'rw24071169processed'} ...
);

for i = 1:height(behav)
    observation_id = behav.ObservationId{i};
    behavior_name = matlab.lang.makeValidName(behav.Behavior{i});

    if isKey(observation_mapping, observation_id)
        whale_obj_name = observation_mapping(observation_id);
        whale_obj = evalin('base', whale_obj_name);

        tstart = behav.tstart(i);
        tstop = behav.tstop(i);

        % Handle behaviors spanning two videos
        if tstop < tstart
            next_seq = behav.seq(i) + 1;
            next_video_idx = find(behav.seq == next_seq & strcmp(behav.ObservationId, observation_id), 1, 'first');
            if ~isempty(next_video_idx)
                tstop = behav.vidst(next_video_idx);
            else
                fprintf('Warning: No next video found for %s, seq %d.\n', observation_id, behav.seq(i));
                continue;
            end
        end

        time_col = whale_obj.Datenum;
        if isnumeric(time_col)
            tstart = datenum(tstart);
            tstop = datenum(tstop);
        elseif isnumeric(tstart)
            tstart = datetime(tstart, 'ConvertFrom', 'datenum');
            tstop = datetime(tstop, 'ConvertFrom', 'datenum');
        end

        start_idx = find(time_col >= tstart, 1, 'first');
        stop_idx = find(time_col <= tstop, 1, 'last');
        stop_idx = stop_idx + 2;

        if ~isempty(start_idx) && ~isempty(stop_idx)
            whale_obj.(behavior_name)(start_idx:stop_idx) = 1;
            assignin('base', whale_obj_name, whale_obj);
            fprintf('Matched: %s with Behavior: %s, StartIdx: %d, StopIdx: %d\n', ...
                behav.MediaFileName{i}, behavior_name, start_idx, stop_idx);
        else
            fprintf('Warning: No matching start or stop index found for %s\n', behav.MediaFileName{i});
        end
    else
        fprintf('Warning: ObservationId %s not found in mapping.\n', observation_id);
    end
end


%clear("behav", "behavior_name", "i", "observation_id", "observation_mapping", ...
    %"start_idx", "stop_idx", "time_col", "tstop", "tstart", "unique_behaviors", ...
    %"whale_obj", "whale_obj_name");
disp("Section 6 - Behaviours Added");


%%
%Section 7 - Logic Check 
%As with the section above this plots one of the PRHs with marked
%surfacings so that you can ensure stitching was done correctly
% Assuming you are plotting from the whale object (e.g., 'rw24071042processed')
whale_obj = mn3009202542processed;  % Replace with the correct object if necessary

% Extract the necessary data
datenum_values = whale_obj.Datenum;  % Time data
depth_values = whale_obj.Depth;      % Depth data
surface_behavior = whale_obj.FirstSurface;  % Replace 'Surface' with the actual column name for surface behavior

% Plot the Datenum vs Depth, inverting the y-axis
figure;
plot(datenum_values, depth_values, 'b-', 'LineWidth', 1.5);  % Plot depth vs time
hold on;
set(gca, 'YDir', 'reverse');  % Invert the y-axis for depth

% Plot points at y = 0 where the first Surface behavior occurs
surface_idx = find(surface_behavior == 1);  % Find indices where 'Surface' == 1
scatter(datenum_values(surface_idx), zeros(size(surface_idx)), 'r', 'filled');  % Plot red points at y = 0

% Customize plot appearance
xlabel('Datenum');
ylabel('Depth (inverted)');
title('Depth vs Time with Surface Behavior at y = 0');
datetick('x', 'keeplimits');  % Format x-axis as datetime if necessary
grid on;
legend('Depth', 'Surface Behavior');

hold off;
clear("datenum_values", "depth_values", "surface_behavior", "surface_idx", ...
    "whale_obj");
disp("Section 7 - Logic Check Complete");

%%
%Data export to csv for visualisation in R
writetable(mn3009202542processed, 'mn3009202542processed.csv');
%%
%Section 8 - Remove behaviour that is not Mouth Open, obviously if you are
%replicating this code to examine another behaviour then change this
%section. The easiest way to generate the below code is to manually open
%each of your PRH's and manually delete the unecessary columns, then you
%can copy and paste the generated code into this section
rw24071042processed = removevars(rw24071042processed, ["LastSurface","Conspecific"]);
rw24071042processed = removevars(rw24071042processed, ["Benthos","Roll","FirstSurface"]);
rw24071344processed = removevars(rw24071344processed, ["Zooplankton","FirstSurface","LastSurface","Conspecific"]);
rw24070745processed = removevars(rw24070745processed, ["Benthos","FirstSurface","LastSurface","Conspecific"]);
rw24071169processed = removevars(rw24071169processed, ["Conspecific","FirstSurface","LastSurface","Roll"]);
disp("Section 8, Unecessary Behaviour Removal Complete"); 
%%
%Section 9 - Combine Data
% Combine the datasets (adjust as needed based on your data structure)
combinedData = [];

for i = 1:length(whale_names)
    whale_name = whale_names(i).name;
    if exist(whale_name, 'var')
        % Get the data associated with the variable name in whale_name
        currentData = eval(whale_name); 
        % Combine the data (assuming vertically stacking is desired)
        combinedData = [combinedData; currentData];
    else
        warning('Variable %s does not exist in the workspace.', whale_name);
    end
end

save('combined.mat', 'combinedData') ;
clear("currentData", "i", "whale_name");
disp("Section 9, Data Combination - Complete");

%%
%Section 10 - Normalisation 
% This section normalises certainfor cross tag comparison
%Similar to other scripts in this package we will normalise certain metrics
%to enable inter-tag comparison 

normData = combinedData;
normData = removevars(normData, ["binary_peaks","DBA","VDBA"]);

%Index the columns to normalise
exclude_columns = find(ismember(normData.Properties.VariableNames, {'Depth', 'Daylight', 'Hour', 'Datenum', 'MouthOpen'}));



% Loop through each column except 'datenum' and 'fluking_event' to normalize
for i = 1:width(normData)
    if ~ismember(i, exclude_columns)
        % Normalize the column by subtracting mean and dividing by standard deviation
        normData{:, i} = (normData{:, i} - mean(normData{:, i}, 'omitna')) / std(normData{:, i}, 'omitna');
    end
end
save('normData.mat', 'normData');
clear("i", "exclude_columns");
disp("Section 10, Normalisation - Complete");


%%
%This section pulls a piece of one of the TDR's which I have previously
%determined based on the video audit. The purpose of this is to build
%figures which will incorporate screenshots from the video audit. 
% Extract the relevant data
figure_data = rw24070745processed(5984:11074, :);
fig_depth = figure_data.Depth;     % Depth values
fig_time = figure_data.Datenum;    % Time values
fig_open = figure_data.MouthOpen;  % MouthOpen status (binary)

% Initialize the figure
figure()
hold on;

% Loop through data points and plot line segments with different colors
for i = 1:length(fig_time) - 1
    % Determine the color based on the MouthOpen status
    if fig_open(i) == 1
        line_color = '#4E8056'; % Fern Greenfor MouthOpen = 1
    else
        line_color = '#D2E4C4'; % Tea Green for MouthOpen = 0
    end
    
    % Plot the line segment
    plot(fig_time(i:i+1), fig_depth(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
ylim([ 0 50])
xlabel('Time (ADT)');
set(gca, 'FontSize', 20);
% Reverse the y-axis to show depth correctly
set(gca, 'YDir', 'reverse');
exportgraphics(gcf, 'tdr_signal.png', 'BackgroundColor', 'none');
hold off;
%%
fig_fluke = figure_data.fluking_signal;
hold on;

% Loop through data points and plot line segments with different colors
for i = 1:length(fig_time) - 1
    % Determine the color based on the MouthOpen status
    if fig_open(i) == 1
        line_color = '#633365'; % Finn for MouthOpen = 1
    else
        line_color = '#DBBADD'; % Pink Lavender for MouthOpen = 0
    end
    
    % Plot the line segment
    plot(fig_time(i:i+1), fig_fluke(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
hold off;

%%
fig_speed = figure_data.speed;
hold on;

% Loop through data points and plot line segments with different colors
for i = 1:length(fig_time) - 1
    % Determine the color based on the MouthOpen status
    if fig_open(i) == 1
        line_color = '#06559F'; % Lapis Lazuli for MouthOpen = 1
    else
        line_color = '#B1CEDF'; % Colombia Blue for MouthOpen = 0
    end
    
    % Plot the line segment
    plot(fig_time(i:i+1), fig_speed(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
hold off;

%%
% Extract relevant data
figure_data = rw24070745processed(5984:11074, :);
fig_depth = figure_data.Depth;     % Depth values
fig_time = figure_data.Datenum;    % Time values
fig_open = figure_data.MouthOpen;  % MouthOpen status (binary)
fig_fluke = figure_data.fluking_signal; % Fluking signal
fig_speed = figure_data.speed;     % Speed

% Initialize figure
figure()

% Define layout with 4 rows to allow depth plot to be larger
tiledlayout(5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Fluking signal plot (top)
ax1 = nexttile(1);
hold on;
for i = 1:length(fig_time) - 1
    % Determine line color based on MouthOpen status
    if fig_open(i) == 1
        line_color = '#633365'; % Finn
    else
        line_color = '#DBBADD'; % Pink Lavender
    end
    plot(fig_time(i:i+1), fig_fluke(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
%ylabel('Fluking Signal');
set(gca, 'xtick', [], 'xticklabel', []);
set(gca, 'FontSize', 20); % Adjust 14 to desired size
hold off;

% Speed plot (middle)
ax2 = nexttile(2);
hold on;
for i = 1:length(fig_time) - 1
    % Determine line color based on MouthOpen status
    if fig_open(i) == 1
        line_color = '#06559F'; % Lapis Lazuli
    else
        line_color = '#B1CEDF'; % Columbia Blue
    end
    plot(fig_time(i:i+1), fig_speed(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
%ylabel('Speed (m/s)');
set(gca, 'xtick', [], 'xticklabel', []);
set(gca, 'FontSize', 20); % Adjust 14 to desired size
ylim([0 10]);
hold off;

% Depth plot (bottom, larger, spanning 2 rows)
ax3 = nexttile(3, [3, 1]); % Depth plot spans two rows
hold on;
for i = 1:length(fig_time) - 1
    % Determine line color based on MouthOpen status
    if fig_open(i) == 1
        line_color = '#4E8056'; % Fern Green
    else
        line_color = '#D2E4C4'; % Tea Green
    end
    plot(fig_time(i:i+1), fig_depth(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
set(ax3, 'YDir', 'reverse'); % Reverse y-axis
ylim([0 50]);
%ylabel('Depth (m)');
xlabel('Time (ADT)');
set(gca, 'FontSize', 20); % Adjust 14 to desired size
hold off;

% Link x-axes
linkaxes([ax1, ax2, ax3], 'x');

% Save the figure
exportgraphics(gcf, 'tdr_three_panel.png', 'BackgroundColor', 'none');

%%
% Separate numeric columns
numericData = combinedData(:, varfun(@isnumeric, combinedData, 'OutputFormat', 'uniform'));


% Extract the MouthOpen column and validate its presence
if ismember('MouthOpen', combinedData.Properties.VariableNames)
    MouthOpen = combinedData.MouthOpen; % Extract the MouthOpen column
else
    error('MouthOpen column not found in combinedData.');
end

% Remove the MouthOpen column from the numeric data
numericData = removevars(numericData, 'MouthOpen');


% Split the numeric data based on MouthOpen values
dataMouthOpen1 = numericData{MouthOpen == 1, :}; % Rows where MouthOpen is 1
dataMouthOpen0 = numericData{MouthOpen == 0, :}; % Rows where MouthOpen is 0

% Now calculate the mean and standard deviation for the transformed data
meanMouthOpen1 = mean(dataMouthOpen1(:,:), 'omitnan');
meanMouthOpen0 = mean(dataMouthOpen0(:,:), 'omitnan');
stdMouthOpen1 = std(dataMouthOpen1(:,:), 'omitnan');
stdMouthOpen0 = std(dataMouthOpen0(:,:), 'omitnan');

% Combine results into tables for clarity
averageResults = table(meanMouthOpen1', meanMouthOpen0', ...
    'VariableNames', {'MouthOpen1', 'MouthOpen0'}, ...
    'RowNames', numericData.Properties.VariableNames);

stdResults = table(stdMouthOpen1', stdMouthOpen0', ...
    'VariableNames', {'MouthOpen1', 'MouthOpen0'}, ...
    'RowNames', numericData.Properties.VariableNames);

%%
%Section 6 - Mouth Open Characteristics 

%This section first creates a large data object similar to the one used in
%the Dive_pattern.mat script, then creates boxplots comparing mouth open
%and closed for each of the predictors based on the the predicted value of
%MouthOpen

all_data = combinedData;


%%
% Add MouthOpen column to nn_data and convert to binary (1 for 'Open', 0 for 'Closed')
nn_data = combinedData;


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
    depthFilter = nn_data.Depth < 5;

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
% Extract data
figure_data = rw24070745processed(5984:11074, :);
fig_time = figure_data.Datenum;    % Time values
fig_depth = figure_data.Depth;     % Depth values
fig_open = figure_data.MouthOpen;  % MouthOpen status (binary)
fig_fluke = figure_data.fluking_signal;
fig_speed = figure_data.speed;

% Initialize figure
figure();

%% Plot Fluking Signal (Top, 1/4 height)
subplot(5,1,1);
hold on;
colors = {'#DBBADD', '#633365'}; % {Pink Lavender, Finn}
for i = 1:length(fig_time) - 1
    line_color = colors{fig_open(i) + 1}; % Select color based on MouthOpen status
    plot(fig_time(i:i+1), fig_fluke(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
ylabel('Fluking Signal (rad)');
set(gca, 'XTickLabel', []); % Remove x-axis labels for clarity
set(gca, 'XTick', []); % Remove x-axis ticks
set(gca, 'FontSize', 12); % Adjust 14 to desired size

hold off;

%% Plot Speed (Middle, 1/4 height)
subplot(5,1,2);
hold on;
colors = {'#B1CEDF', '#06559F'}; % {Colombia Blue, Lapis Lazuli}
for i = 1:length(fig_time) - 1
    line_color = colors{fig_open(i) + 1}; % Select color based on MouthOpen status
    plot(fig_time(i:i+1), fig_speed(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
ylabel('Speed (m/s)');
set(gca, 'XTickLabel', []); % Remove x-axis labels for clarity
set(gca, 'XTick', []); % Remove x-axis ticks
set(gca, 'FontSize', 12); % Adjust 14 to desired size

hold off;

%% Plot Depth (Bottom, 1/2 height)
subplot(5,1,[3 5]); % Make depth take half of the total figure height
hold on;
colors = {'#D2E4C4', '#4E8056'}; % {Tea Green, Fern Green}
for i = 1:length(fig_time) - 1
    line_color = colors{fig_open(i) + 1}; % Select color based on MouthOpen status
    plot(fig_time(i:i+1), fig_depth(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
ylim([0 45]);
set(gca, 'YDir', 'reverse'); % Reverse y-axis for depth
ylabel('Depth (m)');
xlabel('Time (ADT)');
set(gca, 'FontSize', 12); % Adjust 14 to desired size

hold off;

%%
%Crazy detail figure
% Extract data
figure_data = rw24070745processed(5984:11074, :);
fig_time = figure_data.Datenum;    % Time values

% Extract variables for plotting
fig_fluke = figure_data.fluking_signal;
fig_odba = figure_data.ODBA;
fig_vdba = figure_data.VDBA;
fig_jerk = figure_data.jerk;

fig_speed = figure_data.speed;
fig_body_orientation = figure_data.body_orientation;
fig_head = figure_data.head;
fig_roll = figure_data.roll;
fig_fluking_rate = figure_data.fluking_event;

fig_gyr_x = figure_data.Gyro_x;
fig_gyr_y = figure_data.Gyro_y;
fig_gyr_z = figure_data.Gyro_z;
fig_acc_x = figure_data.Acc_x;
fig_acc_y = figure_data.Acc_y;
fig_acc_z = figure_data.Acc_z;
fig_depth = figure_data.Depth;

% Initialize figure
figure();

%% Plot Fluking Signal, ODBA, VDBA, and Jerk (Top, 1/4 height)
subplot(8,1,1);
hold on;
plot(fig_time, fig_fluke, '-', 'Color', '#d55e00', 'LineWidth', 1.5); % Pink Lavender
plot(fig_time, fig_odba, '-', 'Color', '#0072b2', 'LineWidth', 1.5); % Finn
plot(fig_time, fig_vdba, '-', 'Color', '#009e73', 'LineWidth', 1.5); % Colombia Blue

%ylabel('Fluking, ODBA, VDBA');
%legend({'Fluking Signal', 'ODBA', 'VDBA'}, 'Location', 'eastoutside');
set(gca, 'XTickLabel', [], 'XTick', [], 'FontSize', 12);
hold off;
%%
subplot(8,1,2);
hold on;
plot(fig_time, fig_jerk, '-', 'Color', '#d55e00', 'LineWidth', 1.5); % Lapis Lazuli
plot(fig_time, fig_fluking_rate/10, '-', 'Color', '#0072b2', 'LineWidth', 1.5); % Pink Lavender


%ylabel('Jerk, Fluking Rate (/10)');
%legend({'Jerk', 'Fluke Rate/10'}, 'Location', 'eastoutside');
set(gca, 'XTickLabel', [], 'XTick', [], 'FontSize', 12);
hold off;
%% Plot Speed, Body Orientation, Head, and Roll (Middle, 1/4 height)
subplot(8,1,3);
hold on;
plot(fig_time, fig_speed, '-', 'Color', '#d55e00', 'LineWidth', 1.5); % Tea Green
plot(fig_time, fig_body_orientation, '-', 'Color', '#0072b2', 'LineWidth', 1.5); % Fern Green
plot(fig_time, fig_head, '-', 'Color', '#009e73', 'LineWidth', 1.5); % Gold
plot(fig_time, fig_roll, '-', 'Color', '#f0e442', 'LineWidth', 1.5); % Orange Red
%ylabel('Speed, Orientation, Head, Roll');
%legend({'Speed', 'Body Orientation', 'Head', 'Roll'}, 'Location', 'none');
set(gca, 'XTickLabel', [], 'XTick', [], 'FontSize', 10);
hold off;

%% Plot Speed, Body Orientation, Head, and Roll (Middle, 1/4 height)
subplot(8,1,4);
hold on;
plot(fig_time, fig_gyr_x, '-', 'Color', '#d55e00', 'LineWidth', 1.5); % Tea Green
plot(fig_time, fig_gyr_y, '-', 'Color', '#0072b2', 'LineWidth', 1.5); % Fern Green
plot(fig_time, fig_gyr_z, '-', 'Color', '#009e73', 'LineWidth', 1.5); % Gold
%ylabel('Gyroscopes');
%legend({'X', 'Y' , 'Z'}, 'Location', 'eastoutside');
set(gca, 'XTickLabel', [], 'XTick', [], 'FontSize', 10);
hold off;
%% Plot Speed, Body Orientation, Head, and Roll (Middle, 1/4 height)
subplot(8,1,5);
hold on;
plot(fig_time, fig_acc_x, '-', 'Color', '#d55e00', 'LineWidth', 1.5); % Tea Green
plot(fig_time, fig_acc_y, '-', 'Color', '#0072b2', 'LineWidth', 1.5); % Fern Green
plot(fig_time, fig_acc_z, '-', 'Color', '#009e73', 'LineWidth', 1.5); % Gold
%ylabel('Accelerometer');
%legend({'X', 'Y' , 'Z'}, 'Location', 'eastoutside');
set(gca, 'XTickLabel', [], 'XTick', [], 'FontSize', 10);
hold off;
%% Plot Depth (Bottom, 1/2 height)
subplot(8,1,[6 8]); % Make depth take half of the total figure height
hold on;
plot(fig_time, fig_depth, '-', 'Color', '#d55e00', 'LineWidth', 1.5); % Fern Green
ylim([0 45]);
set(gca, 'YDir', 'reverse'); % Reverse y-axis for depth
%ylabel('Depth (m)');
%xlabel('Time (ADT)');
set(gca, 'FontSize', 12);
hold off;
