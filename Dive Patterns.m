%% Dive Patterns
%This script visualises dive depth through a 24h cycle
%This requires the all_data output from the Dive_visualise.m script 

data_backup = load('all_data.mat');

data_backup = data_backup.all_data{:, {'Datenum', 'Depth'}};

% Step 2: Filter the data so that Depth does not go above 60
filtered_data = array2table(data_backup);

filtered_data.Date = datetime(filtered_data.(filtered_data.Properties.VariableNames{1}), 'ConvertFrom', 'datenum');
filtered_data.Properties.VariableNames(1) = "Date";
filtered_data.Properties.VariableNames(2) = "Depth";

% Now you can use the hour() function on the converted datetime
filtered_data.Hour = hour(filtered_data.Date);



% Step 4: Define bin edges for hour and depth
% Hour binning remains 0-23, Depth binned between 0 and 60
hour_bins = 0:23;  % 24 hours of the day
depth_bins = linspace(0, 60, 30);  % 100 depth bins between 0 and 60 meters

% Step 5: Bin the data and calculate density
% Use histcounts2 to count the number of observations in each bin
[counts, hour_edges, depth_edges] = histcounts2(filtered_data.Hour, filtered_data.Depth, hour_bins, depth_bins);

% Step 6: Create the heatmap
figure;
imagesc(hour_edges(1:end-1), depth_edges(1:end-1), log(counts)');  % Transpose counts to align with axes
% Step 7: Customize the plot
colorbar;  % Show colorbar for density
xlabel('Hour of Day');
ylabel('Depth (0-60)');
title('Density of Observations by Hour and Depth (Up to 60m)');
colormap('parula');  % Use a color map for better visualization