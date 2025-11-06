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

    % Call the getgape function (without providing lnth)
    [gapearea, ci_gapes, width, ci_width, lnth] = getgape(age);

    % Store the results in the corresponding columns
    WhaleAge.gapearea(i) = gapearea;
    WhaleAge.ci_gapes(i) = ci_gapes;
    WhaleAge.width(i) = width;
    WhaleAge.ci_width_lower(i) = ci_width(1); % Store lower bound of ci_width
    WhaleAge.ci_width_upper(i) = ci_width(2); % Store upper bound of ci_width
    WhaleAge.lnth(i) = lnth;
end

save('WhaleAge.mat', "WhaleAge");
