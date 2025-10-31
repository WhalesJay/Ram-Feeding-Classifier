%%Visualisation script
%This script sets up the interactive filtering UI for initial review of
%dives. This script uses outputs from the Dive_analyse.m script which needs
%to be run prior to this. 
%In section 3 there is a switch for normalising data, which will normalise
%speed, fluking frequency, jerk and DBA. If you would prefer to plot
%without this please change the switch to false

%%
%Section 1: Remove variables as needed and combine the data into a single
%object.
% Part 1: Load and Store whale_data Structures
whale_names = dir('*processed.mat');
name_list = {}; % Initialize a cell array to store file names

% Loop through the list and filter out files with underscores in their names
for i = 1:length(whale_names)
    % Check if the file name does not contain an underscore
    if isempty(strfind(whale_names(i).name, '_'))
        % Load the file if it meets the criteria
        data = load(whale_names(i).name);
        
        % Create a variable name without the .mat extension
        var_name = erase(whale_names(i).name, '.mat');
        
        % Store the name in name_list for later processing
        name_list{end + 1} = var_name;
        
        % Assign the whale_data structure to the base workspace
        assignin('base', var_name, data.whale_data);
    end
end

% Optionally save name_list to the workspace for later use
clear("whale_names", "var_name", "i", "data");
assignin('base', 'whale_names', name_list);
clear("name_list");

% Part 2: Concatenate Data from All whale_data Objects
columns_to_use = {'Datenum', 'binary_peaks', 'jerk', 'Depth', 'DBA', 'speed', 'fluking_event', 'head', 'pitch', 'roll'};
all_data = table(); % Initialize an empty table to hold data from all whale names

% Loop through each whale_name object and concatenate the specified columns
for i = 1:length(whale_names)
    whale_name = whale_names{i}; % Get the whale name without .mat
    whale_data = evalin('base', whale_name); % Access the whale data object from the base workspace
    
    % Check if all required columns exist
    if all(ismember(columns_to_use, whale_data.Properties.VariableNames))
        current_data = whale_data(:, columns_to_use); % Extract relevant data as a table
        
        % Add a new column for the whale name
        current_data.WhaleName = repmat(string(whale_name), height(current_data), 1);
        all_data = [all_data; current_data]; % Append to all_data table
    else
        fprintf('Warning: %s does not contain all required columns.\n', whale_name);
    end
end


% Remove rows with NaN values in any column
all_data = rmmissing(all_data); % This will remove any rows containing NaNs
save('all_data.mat', 'all_data');
clear('data', 'i', 'var_name');
disp('Section 1 Complete - All Data Created');
%%
%Section 2: Normalising for cross tag comparison, skip this section
% Assuming all_data is your data table or matrix
% Find the indices of the 'datenum' and 'fluking_event' columns
exclude_columns = find(ismember(all_data.Properties.VariableNames, {'Datenum', 'fluking_event', 'WhaleName', 'Depth'}));

% Initialize norm_data as a copy of all_data
norm_data = all_data;

% Loop through each column except 'datenum' and 'fluking_event' to normalize
for i = 1:width(all_data)
    if ~ismember(i, exclude_columns)
        % Normalize the column by subtracting mean and dividing by standard deviation
        norm_data{:, i} = (all_data{:, i} - mean(all_data{:, i})) / std(all_data{:, i});
    end
end

norm_data.Datenum = datetime(all_data.Datenum, 'ConvertFrom', 'datenum');
disp('Section 2 - Normalisation Complete')
%%
%This is Jays filtering plot function on #MEGADATAFRAME
%For some reason the axis invert is not working, I don't know why, but
%personally I don't care, if you do hang upside down like batman and
%interpret from there. If you would like to plot the non-normalised data
%then change the line below to all_data rather than norm_data
plotWhaleDataWithUI(norm_data)
