


%% Section 1: Load Your Data
% Load your dataset here
% Assuming 'features' is a table or array containing the input data
% Assuming 'labels' is a table or array containing the output data (MouthOpen)
% For example, you can load your data like this:
% load('yourDataset.mat'); % Adjust this line based on your actual dataset

%% Section 2: Split Data into Training and Testing Sets
% Split the data into training and testing sets
cv = cvpartition(size(features, 1), 'HoldOut', 0.3); % 70% train, 30% test
X_train = features(training(cv), :);
Y_train = labels(training(cv), :);
X_test = features(test(cv), :);
Y_test = labels(test(cv), :);

%% Section 3: Prepare Y_train
% If Y_train is a table and 'MouthOpen' is the column name
Y_train = Y_train.MouthOpen; % Extract the 'MouthOpen' column

% Convert to categorical (if not already categorical)
Y_train = categorical(Y_train);

% Calculate the number of features and classes
numFeatures = size(X_train,2);
numClasses = numel(categories(Y_train));


%% Section 4: Prepare X_train for LSTM
% Extract features into an array
X_data = table2array(X_train);  % Convert features to an array

% Get the total number of rows in X_data
numRows = size(X_data, 1);  % Total number of rows in the dataset

% Define numTimesteps and calculate the number of sequences
numTimesteps = 100;  % Set the number of timesteps per sequence
numSequences = floor(numRows / numTimesteps);  % Calculate number of complete sequences

% Remove excess rows if they don't fit into complete sequences
X_data = X_data(1:(numSequences * numTimesteps), :);  % Ensure X_data fits complete sequences

% Reshape the data: [numTimesteps, numFeatures, numSequences]
X_data = reshape(X_data, [numTimesteps, size(X_data, 2), numSequences]);  % Reshape into 3D

% Convert into N-by-1 cell array, where each cell contains a T-by-F sequence
X_train = squeeze(mat2cell(X_data, numTimesteps, size(X_data, 2), ones(1, numSequences)));  % Convert to cell array

%% Section 5: Prepare Y_train for LSTM
% Ensure Y_data is structured correctly
Y_data = Y_train(:);  % Convert to a column vector if necessary

% Remove excess rows if they don't fit into complete sequences
Y_data = Y_data(1:(numSequences * numTimesteps), :);  

% Reshape into N sequences of T timesteps with F features (F = 1 for MouthOpen)
Y_data = reshape(Y_data, [numTimesteps, 1, numSequences]);  % Here F = 1

% Convert into N-by-1 cell array, where each cell contains a T-by-F sequence
Y_train = squeeze(mat2cell(Y_data, numTimesteps, 1, ones(1, numSequences))); 

% Check the final sizes of X_train and Y_train
disp(['Size of X_train: ', num2str(size(X_train))]);  % Should show [numSequences, 1]
disp(['Size of Y_train: ', num2str(size(Y_train))]);  % Should show [numSequences, 1]

%% Section 6: Define LSTM Layers
% Define LSTM layers and train
inputSize = 15;  % Number of features
numHiddenUnits = 100;

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits, 'OutputMode', 'last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

%% Section 7: Set Training Options
options = trainingOptions('adam', ...
    'MaxEpochs', 30, ...
    'MiniBatchSize', 64, ...
    'ValidationData', {X_test, Y_test}, ...
    'Plots', 'training-progress');

%% Section 8: Train the Network
% Train the network
net = trainnet(X_train, Y_train, layers, options);
