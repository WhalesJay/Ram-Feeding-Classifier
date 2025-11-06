
%Load in the tagged whale information
WhaleAge = readtable("C:\Users\Jay\OneDrive - Dalhousie University\Publications\Feeding Kinematics\Gape\Whale_Age.csv")


% Outputs from getgape will be saved as new columns in WhaleAge

% Initialize new columns in the table
WhaleAge.gapearea = NaN(height(WhaleAge), 1);
WhaleAge.ci_gapes = NaN(height(WhaleAge), 1);
WhaleAge.width = NaN(height(WhaleAge), 1);
WhaleAge.ci_width_lower = NaN(height(WhaleAge), 1); % Lower bound of ci_width
WhaleAge.ci_width_upper = NaN(height(WhaleAge), 1); % Upper bound of ci_width
WhaleAge.lnth = NaN(height(WhaleAge), 1);

% Iterate through each row of the table and calculate values
for i = 1:height(WhaleAge)
    % Extract the age value for the current row
    age = WhaleAge.Age(i);

    % If age is greater than 30, set it to 30
    if age > 30
        age = 30;
    end

    % Call the getgape function (without providing lnth)
    [gapearea, ci_gapes, width, ci_width, lnth] = getgape(age);

    % Store the results in the corresponding columns
    %I'm a dope and I calculated area as m*cm so the area comes out way too
    %big for m^2. I know how, I just don't know where in the code to
    %correct because this is partially not my code, so below is an odd
    %looking sticking plaster

    WhaleAge.gapearea(i) = gapearea/100;
    WhaleAge.ci_gapes(i) = ci_gapes/100;
    WhaleAge.width(i) = width;
    WhaleAge.ci_width_lower(i) = ci_width(1); % Store lower bound of ci_width
    WhaleAge.ci_width_upper(i) = ci_width(2); % Store upper bound of ci_width
    WhaleAge.lnth(i) = lnth;
end


save('WhaleAge.mat', "WhaleAge");
%%
% Assuming the following:
% - whale_names is a cell array of object names (strings)
% - Each object in whale_names contains 'MouthOpen' and 'speed' columns
% - WhaleAge table has columns 'ID', 'gapearea', and 'volume_filtered'
%load("WhaleAge.mat")
% Initialize the 'volume_filtered' column in WhaleAge
WhaleAge.volume_filtered = zeros(height(WhaleAge), 1);

% Loop through each whale name
for i = 1:length(whale_names)
    % Get the current whale name
    whale_name = whale_names(i).name;
    
    whale_data = eval(whale_name); 
    speed = whale_data.nnspeed;
    MouthOpen = whale_data.MouthOpen;
    % Calculate the total distance traveled when the mouth is open
    % - nrow(MouthOpen == 1) gives the total number of rows with MouthOpen = 1
    % - Average speed is calculated only during mouth open
    mouth_open_indices = MouthOpen == 1;
    
    avg_speed = mean(speed(mouth_open_indices), 'omitnan'); % Average speed when mouth is open
    nrows_open = sum(MouthOpen); % Number of rows with MouthOpen = 1
    abs_open = nrows_open /10; % Convert rows to seconds (10 Hz)
    % Calculate the total distance traveled (in meters)
    distance_traveled = avg_speed * abs_open; %Average speed (m/s) * time open (s) = distance (m)
    
    % Get the gape area for the current whale from WhaleAge
    % Match the whale ID in WhaleAge to the current whale name
    whale_row = find(strcmp(WhaleAge.ID, whale_name));
    if isempty(whale_row)
        warning('Whale ID %s not found in WhaleAge. Skipping.', whale_name);
        continue;
    end
    gape_area = WhaleAge.gapearea(whale_row);

    % Calculate the volume filtered and update WhaleAge
    volume_filtered = gape_area * distance_traveled; % Volume in m³ (or appropriate units)
    WhaleAge.volume_filtered(whale_row) = (volume_filtered);
    WhaleAge.Abs_time(whale_row) = abs_open;
    % Calculate the confidence interval for the filtered volume
    ci_filtered_volume = WhaleAge.ci_gapes(whale_row) * distance_traveled;
    
    % Add the confidence interval columns to WhaleAge
    WhaleAge.ci_volume_filtered(whale_row) = ci_filtered_volume;
    
    % Convert TagDuration to hours
    tag_duration_hours = hours(WhaleAge.TagDuration(whale_row));
    WhaleAge.avg_speed(whale_row) = avg_speed;
    WhaleAge.dist(whale_row) = distance_traveled ; 
    % Calculate volume filtered per hour
    WhaleAge.volume_filtered_per_hour(whale_row) = WhaleAge.volume_filtered(whale_row) / tag_duration_hours;
    percentage_open = (sum(mouth_open_indices) / length(MouthOpen)) * 100; % Percentage
    % Add the percentage to WhaleAge
    WhaleAge.percentage_mouth_open(whale_row) = percentage_open;
end

%%
% Create a figure
% Sort the WhaleAge table by the percentage of time feeding (from smallest to largest)
[~, sortedIdx] = sort(WhaleAge.percentage_mouth_open);

WhaleAge.Age = double(WhaleAge.Age); % For text-to-number conversion


% Reorder WhaleAge according to the sorted index
WhaleAge = WhaleAge(sortedIdx, :);
figure;

% Create the bar plot for Volume Filtered per Hour
b = bar(WhaleAge.volume_filtered_per_hour);
hold on;

% Normalize the percentage time with mouth open to a [0, 1] range for coloring
normalized_percentage = WhaleAge.percentage_mouth_open / 100; 

% Apply the normalized percentage to the bar color fill
b.FaceColor = 'flat';  % Enable individual bar coloring
b.CData = [normalized_percentage, zeros(height(WhaleAge), 1), ones(height(WhaleAge), 1)]; % RGB color map

% Add labels and title
xlabel('Name / Catalogue ID');
ylabel('Volume Filtered per Hour (m^3/hour)');


% Set x-ticks to each whale and rotate labels
xticks(1:height(WhaleAge)); 
xticklabels(WhaleAge.IndividualID_Name_Catlogue__); % Set the x-tick labels to whale names or IDs
xtickangle(45); % Rotate labels for better readability

% Add a colorbar for the percentage color mapping
colorbar;
colormap('cool'); % Change colormap to 'cool' (or choose any other like 'jet', 'parula', etc.)
caxis([0 1]); % Color axis limits for percentage range (0 to 1)

% Add textboxes with the percentage feeding value on each bar
for i = 1:height(WhaleAge)
    % Display percentage of time with mouth open as text on each bar
    text(i, WhaleAge.volume_filtered_per_hour(i), ...
        sprintf('%.1f%%', WhaleAge.percentage_mouth_open(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');   
end
% Add textboxes with the percentage feeding value on each bar
for i = 1:height(WhaleAge)
    % Display age as text inside the bottom of each bar
    text(i, WhaleAge.Age(i) * 0.05, num2str(WhaleAge.Age(i)), ... % Position slightly above the bottom
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');   
end
ax = gca; % Get current axes
ax.TickLength = [0 0]; % Set tick length to zero
% Adjust plot appearance
hold off;

%%

% Sort the WhaleAge table by the percentage of time feeding (from smallest to largest)
[~, sortedIdx] = sort(WhaleAge.percentage_mouth_open);

% Reorder WhaleAge according to the sorted index
WhaleAge_sorted = WhaleAge(sortedIdx, :);

WhaleAge_sorted.EG = WhaleAge_sorted.IndividualID_Name_Catlogue__;

%Because the above line is for a trial I am just removing the names
%manually, if I am rewriting this code for a publication or you are
%unfortunate enoough to read this line it means I haven't added EG number
%as a column in the input. Either add it back or do it the lazy (higher
%effort way) like me
% Custom Colors (Change these to your desired hex codes)
feeding_color = '#1f77b4';  % Blue color for feeding
nonfeeding_color = '#d62728';  % Red color for non-feeding

% Create figure
figure;

% Loop through each whale and plot the bars
hold on;

% Set positions for the bars (X-axis), now sorted
bar_positions = 1:height(WhaleAge_sorted); 

% Plot bars showing volume filtered per hour
hBar = bar(bar_positions, WhaleAge_sorted.volume_filtered_per_hour, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.5);
hBar.FaceColor = [0.5, 0.5, 0.5];  % Set the bar color to grey initially

% Loop through each whale to overlay text for percentage feeding
for i = 1:height(WhaleAge_sorted)
    % Get the volume filtered per hour for the current whale
    volume = WhaleAge_sorted.volume_filtered_per_hour(i);
    
    % Get the percentage time feeding (mouth open) for the current whale
    feeding_percentage = WhaleAge_sorted.percentage_mouth_open(i);
    
    % Add the percentage feeding value as text above the bar (without "feeding")
    text(bar_positions(i), volume + 0.05, sprintf('%.1f%%', feeding_percentage), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', 'black');
end

% Set labels and title
xlabel('Individual ID');
ylabel('Volume Filtered per Hour (m^3/hour)');
title('Volume Filtered per Hour with Percentage Time Mouth Open');

% Set x-ticks to each whale and rotate labels for better readability
xticks(1:height(WhaleAge_sorted));
xticklabels(WhaleAge_sorted.EG);  % Set the x-tick labels to whale names or IDs
xtickangle(45);  % Rotate labels for better readability

% Adjust plot appearance
hold off;

%%
% Sort the WhaleAge table by the age (smallest to largest)
[~, sortedIdx] = sort(WhaleAge.Age);  % Sort by age
WhaleAge_sorted = WhaleAge(sortedIdx, :);

% Remove the 11th row from WhaleAge_sorted
WhaleAge_sorted(11, :) = [];

WhaleAge_sorted.abstime_open_hours = WhaleAge_sorted.Abs_time / 3600;

% Create a figure for the scatter plot
figure;
hold on;

% Set the figure and axes background to white
set(gcf, 'Color', 'w');  % Figure background
set(gca, 'Color', 'w');  % Axes background

% Use WhaleAge_sorted.Age as the color data
feeding_time = WhaleAge_sorted.abstime_open_hours;  % x-axis
volume_filtered = WhaleAge_sorted.volume_filtered_per_hour;  % y-axis
point_color = WhaleAge_sorted.Age;  % Color data based on age

% Plot the scatter points
scatter(feeding_time, volume_filtered, 100, point_color, 'filled', 'MarkerEdgeColor', 'k');

% Set colormap and add a colorbar for age
colormap(pink);  
c = colorbar;
c.Ticks = [1, 10, 20, 30, 40, 50];   % Dynamic tick values
c.Label.String = 'Age';

% Add axis labels
xlabel('Absolute Time Mouth Open (hours)');
ylabel('Volume Filtered per Hour (m^3/hour)');

% Adjust plot appearance
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid off;
box on;
hold off;

%%

% Create a figure for the scatter plot
figure;
hold on;

% Set the figure and axes background to white
set(gcf, 'Color', 'w');  % Figure background
set(gca, 'Color', 'w');  % Axes background

% Use WhaleAge_sorted.Age as the color data
feeding_time = WhaleAge_sorted.percentage_mouth_open;  % x-axis
volume_filtered = WhaleAge_sorted.volume_filtered_per_hour;  % y-axis
point_color = WhaleAge_sorted.Age;  % Color data based on age

% Plot the scatter points
scatter(feeding_time, volume_filtered, 100, point_color, 'filled', 'MarkerEdgeColor', 'k');

% Set colormap and add a colorbar for age
colormap("parula");  
c = colorbar;
c.Ticks = [1, 10, 20, 30, 40, 50];   % Dynamic tick values
c.Label.String = 'Age';

% Add axis labels
xlabel('Percentage of time with mouth open');
ylabel('Volume Filtered per Hour (m^3/hour)');

% Adjust plot appearance
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid off;
box on;
hold off;

%%
%WhaleAge_sorted = WhaleAge(sortedIdx, :);
%WhaleAge_sorted(11,:) = []; %Temporary line to remove the unprocessed tag 
%WhaleAge(1,:) = []; %As above
% Create a figure for the scatter plot
figure;
hold on;

% Set the background to white
set(gcf, 'Color', 'w');  % Figure background
set(gca, 'Color', 'w');  % Axes background


% Define a colormap for age (you can customize the colormap as needed)
age_colormap = parula(height(WhaleAge_sorted));  % Parula colormap scaled to the number of whales
normalized_ages = WhaleAge_sorted.Age;
colormap_size = size(age_colormap, 1);
scaled_colormap_indices = round(normalized_ages * (colormap_size - 1)) + 1;  % Offset by 1 for 1-based indexing
age_colormap = age_colormap(scaled_colormap_indices, :);

% Loop through each whale and plot the points
for i = 1:height(WhaleAge_sorted)
    whale_name = WhaleAge_sorted(i, 1);  % Convert the ID to string
    whale_name = whale_name.ID{1};
    whale_name = char(whale_name);
    whale_data = eval(whale_name);  % Use eval to get the corresponding data
    WhaleAge_sorted.abstime_open(i) = (WhaleAge_sorted.percentage_mouth_open(i)/100) * size(whale_data, 1);
    WhaleAge_sorted.abstime_open_hours = WhaleAge_sorted.abstime_open / 3600;
    % Extract data for the current whale
    feeding_time = WhaleAge_sorted.abstime_open(i);
    volume_filtered = WhaleAge_sorted.volume_filtered_per_hour(i);
    age = WhaleAge_sorted.Age(i);  % Assuming 'Age' is a numeric column in WhaleAge_sorted
    
    % Use age to set the color of the point (scaled index)
    point_color = age_colormap(i, :);  % Select color based on row index
    
    % Plot the scatter point
    scatter(feeding_time, volume_filtered, 100, point_color, 'filled', 'MarkerEdgeColor', 'k');
end

% Add a colorbar to indicate age
colormap(age_colormap);
c = colorbar;
c.Label.String = 'Age';

% Add axis labels (no title)
xlabel('Absolute Time Mouth Open');
ylabel('Volume Filtered per Hour (m^3/hour)');

% Adjust plot appearance
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid off;
box on;
hold off;

%%
%Feed percent vs Age

%WhaleAge_sorted = WhaleAge(sortedIdx, :);
%WhaleAge_sorted(11,:) = []; %Temporary line to remove the unprocessed tag 
%WhaleAge(1,:) = []; %As above
% Create a figure for the scatter plot
figure;
hold on;

% Set the background to white
set(gcf, 'Color', 'w');  % Figure background
set(gca, 'Color', 'w');  % Axes background


% Define a colormap for age (you can customize the colormap as needed)
age_colormap = parula(height(WhaleAge_sorted));  % Parula colormap scaled to the number of whales
normalized_ages = WhaleAge_sorted.volume_filtered_per_hour;
colormap_size = size(age_colormap, 1);
scaled_colormap_indices = round(normalized_ages * (colormap_size - 1)) + 1;  % Offset by 1 for 1-based indexing
age_colormap = age_colormap(scaled_colormap_indices, :);

% Loop through each whale and plot the points
for i = 1:height(WhaleAge_sorted)
    whale_name = WhaleAge_sorted(i, 1);  % Convert the ID to string
    whale_name = whale_name.ID{1};
    whale_name = char(whale_name);
    whale_data = eval(whale_name);  % Use eval to get the corresponding data
    WhaleAge_sorted.abstime_open(i) = (WhaleAge_sorted.percentage_mouth_open(i)/100) * size(whale_data, 1);
    WhaleAge_sorted.abstime_open_hours = WhaleAge_sorted.abstime_open / 3600;
    % Extract data for the current whale
    feeding_time = WhaleAge_sorted.abstime_open(i);
    volume_filtered = WhaleAge_sorted.volume_filtered_per_hour(i);
    age = WhaleAge_sorted.Age(i);  % Assuming 'Age' is a numeric column in WhaleAge_sorted
    
    % Use age to set the color of the point (scaled index)
    point_color = age_colormap(i, :);  % Select color based on row index
    
    % Plot the scatter point
    scatter(feeding_time, volume_filtered, 100, point_color, 'filled', 'MarkerEdgeColor', 'k');
end

% Add a colorbar to indicate age
colormap(age_colormap);
c = colorbar;
c.Label.String = 'Age';

% Add axis labels (no title)
xlabel('Absolute Time Mouth Open');
ylabel('Volume Filtered per Hour (m^3/hour)');

% Adjust plot appearance
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid off;
box on;
hold off;
%%
%%
%A simpler version of previous code 
scatter(WhaleAge.percentage_mouth_open, WhaleAge.volume_filtered_per_hour, 200, WhaleAge.Age, 'filled');
colorbar; % Add a color bar to indicate the age scale
xlabel('Percentage time with mouth open');
ylabel('Volume filtered per hour (m^3/hour)');
title('');
colormap jet; % Set a colormap (you can change to other colormaps like 'parula', 'hot', etc.)
grid off;
set(gca, 'FontSize', 25, 'LineWidth', 1.5);
c = colorbar;
c.Label.String = 'Age';

c.Label.FontSize = 30; % Optionally sets the font size of the colorbar title.
c.Label.Rotation = 0;  % Rotate title by 90 degrees
c.Label.Position = [1.5, 54, 0];  % Adjust title position to be next to the colorbar
%colorbar('off');

%%
%A simpler version of previous code Age vs Feed Percent
scatter(WhaleAge.Age, WhaleAge.percentage_mouth_open, 200, WhaleAge.volume_filtered_per_hour, 'filled');
colorbar; % Add a color bar to indicate the age scale
xlabel('Age (Years)');
ylabel('Percentage time with mouth open');
title('');
colormap jet; % Set a colormap (you can change to other colormaps like 'parula', 'hot', etc.)
grid off;
set(gca, 'FontSize', 25, 'LineWidth', 1.5);
c = colorbar;
c.Label.String = 'Volume filtered per hour (m^3/hour)';

c.Label.FontSize = 20; % Optionally sets the font size of the colorbar title.
c.Label.Rotation = 90;  % Rotate title by 90 degrees
c.Label.Position = [5, 3200, 0];  % Adjust title position to be next to the colorbar
%colorbar('off');

%%
%%
% Assuming whale_names and the selected random whale object exist
% Example: whale_names = ["whale_1", "whale_2", ...];
random_index = randi(length(whale_names));
random_whale_name = whale_names(random_index).name;

% Evaluate the random whale object
random_whale_obj = eval(random_whale_name);

% Extract necessary data
datenum_values = random_whale_obj.Datenum;     % Time data
depth_values = random_whale_obj.Depth;        % Depth data
mouth_open = random_whale_obj.MouthOpen;      % MouthOpen (1 or 0)
speed = random_whale_obj.speed;               % Speed data
fluking_signal = random_whale_obj.fluking_signal; % Fluking signal data

% Identify continuous chunks where MouthOpen == 1
open_indices = find(mouth_open == 1); % Indices where MouthOpen == 1
diff_indices = diff(open_indices);   % Difference between consecutive indices

% Find continuous segments
breaks = [0; find(diff_indices > 1); length(open_indices)];
segments = arrayfun(@(x) open_indices(breaks(x)+1:breaks(x+1)), 1:length(breaks)-1, 'UniformOutput', false);

% Filter segments with at least 600 rows
valid_segments = segments(cellfun(@length, segments) >= 1200);

% Ensure there are valid segments
if isempty(valid_segments)
    error("No continuous chunks with at least 600 rows where MouthOpen == 1.");
end

% Randomly select one valid segment
selected_segment = valid_segments{randi(length(valid_segments))};

% Select the first 600 rows from the chosen segment
selected_indices = selected_segment(1:1200);

% Extract the corresponding data
selected_datenum = datenum_values(selected_indices);
selected_depth = depth_values(selected_indices);
selected_speed = speed(selected_indices);
selected_fluking = fluking_signal(selected_indices);

% Plot the data
figure;

% Subplot 1: Fluking Signal
subplot(3, 1, 1);
plot(selected_datenum, selected_fluking, '-o');
xlabel('Time');
ylabel('Fluking Signal');
title('Fluking Signal vs Time');
grid on;

% Subplot 2: Speed
subplot(3, 1, 2);
plot(selected_datenum, selected_speed, '-o');
xlabel('Time');
ylabel('Speed (m/s)');
title('Speed vs Time');
grid on;

% Subplot 3: Depth
subplot(3, 1, 3);
plot(selected_datenum, selected_depth, '-o');
xlabel('Time');
ylabel('Depth (m)');
title('Depth vs Time');
set(gca, 'YDir', 'reverse'); % Invert y-axis for depth
grid on;

disp("Random continuous 600-row chunk plotted for predicted mouth open: " + random_whale_name);
%%
% Extract the rows where MouthOpen is 1
openRows = rw24071344processed.MouthOpen == 1;

% Calculate the average speed for those rows
avgSpeed = mean(rw24071344processed.speed(openRows), 'omitna');

% Count the number of rows where MouthOpen is 1
numRows = sum(openRows);

% Print the results
fprintf('Average Speed when MouthOpen = 1: %.2f\n', avgSpeed);
fprintf('Number of Rows where MouthOpen = 1: %d\n', numRows);

