function Y = dive_stats(P, X, dive_cues, sampling_rate, prop, angular, X_name, na_rm)
    if nargin < 4
        error('Not enough input arguments. Please provide P, X, dive_cues, and sampling_rate.');
    end

    if nargin < 5
        prop = 0.85;  % Default proportion
    end
    
    if nargin < 6
        angular = true;  % Default angular data
    end
    
    if nargin < 7
        X_name = 'pitch';  % Default name for X
    end
    
    if nargin < 8
        na_rm = true;  % Default for na.rm
    end

    % Convert dive_cues from table to numeric array and multiply by sampling rate
    start_times = table2array(dive_cues(:, 1));  % Extract start times from the first column
    end_times = table2array(dive_cues(:, 2));    % Extract end times from the second column
    di = round([start_times, end_times] * sampling_rate);   % Multiply by sampling rate and round the values

    % Initialize output dataframe
    Y = table((1:size(dive_cues, 1))', 'VariableNames', {'num'});
    
    % Loop through each dive/flight
    for d = 1:size(dive_cues, 1)
        z = P(di(d, 1):di(d, 2));  % Get depth data for the current dive
        Y.max(d) = max(z, [], 'omitnan');
        
        % Find the phase based on prop
        pt = find(z > prop * max(z), 1, 'first'):find(z > prop * max(z), 1, 'last');
        Y.dur(d) = dive_cues{d, 2} - dive_cues{d, 1};
        
        % Calculate destination phase times
        Y.dest_st(d) = pt(1) / sampling_rate + dive_cues{d, 1};
        Y.dest_et(d) = pt(end) / sampling_rate + dive_cues{d, 1};
        Y.dest_dur(d) = Y.dest_et(d) - Y.dest_st(d);
        
        % Calculate to and from phases
        Y.to_dur(d) = pt(1) / sampling_rate;
        Y.from_dur(d) = (length(z) - pt(end)) / sampling_rate;
        
        % Calculate vertical rates
        Y.to_rate(d) = (z(pt(1)) - z(1)) / Y.to_dur(d);
        Y.from_rate(d) = (z(end) - z(pt(end))) / Y.from_dur(d);
        
        % Process angular data if provided
        if ~isempty(X)
            if angular  % If angular data (e.g., pitch)
                a = X(di(d, 1):di(d, 2));  % Get angular data for the current dive
                at = a(1:pt(1));  % Data for the "to" phase
                af = a(pt(end):end);  % Data for the "from" phase
                ad = a(pt(1):pt(end));  % Data for the "destination" phase
                
                % Remove NaNs if necessary
                if na_rm
                    a = a(~isnan(a));
                    at = at(~isnan(at));
                    af = af(~isnan(af));
                    ad = ad(~isnan(ad));
                end
                
                % Compute mean and variance for the angles
                Y.mean_angle(d) = mean(a, 'omitnan');
                Y.angle_var(d) = var(a, 'omitnan');
                Y.mean_to_angle(d) = mean(at, 'omitnan');
                Y.mean_dest_angle(d) = mean(ad, 'omitnan');
                Y.mean_from_angle(d) = mean(af, 'omitnan');
                Y.to_angle_var(d) = var(at, 'omitnan');
                Y.dest_angle_var(d) = var(ad, 'omitnan');
                Y.from_angle_var(d) = var(af, 'omitnan');
            else  % If non-angular data
                a = X(di(d, 1):di(d, 2));  % Get the data for the current dive
                at = a(1:pt(1));  % Data for the "to" phase
                af = a(pt(end):end);  % Data for the "from" phase
                ad = a(pt(1):pt(end));  % Data for the "destination" phase
                
                % Compute mean and standard deviation for each phase
                Y.mean_aux(d) = mean(a, 'omitnan');
                Y.aux_sd(d) = std(a, 'omitnan');
                Y.mean_to_aux(d) = mean(at, 'omitnan');
                Y.mean_dest_aux(d) = mean(ad, 'omitnan');
                Y.mean_from_aux(d) = mean(af, 'omitnan');
                Y.to_aux_sd(d) = std(at, 'omitnan');
                Y.dest_aux_sd(d) = std(ad, 'omitnan');
                Y.from_aux_sd(d) = std(af, 'omitnan');
            end
        end
    end
    
    % Rename the columns of Y if needed
    if ~strcmp(X_name, 'pitch')
        Y.Properties.VariableNames = strrep(Y.Properties.VariableNames, 'angle', X_name);
        Y.Properties.VariableNames = strrep(Y.Properties.VariableNames, 'aux', X_name);
    end
    
    % Add start and end times from dive_cues to the output table
    Y.st = dive_cues{:, 1};
    Y.et = dive_cues{:, 2};
    
    % Remove any variables that were not created (due to missing data or conditions)
    Y = removevars(Y, setdiff(Y.Properties.VariableNames, ...
        {'num', 'max', 'st', 'et', 'dur', 'dest_st', 'dest_et', 'dest_dur', ...
         'to_dur', 'from_dur', 'mean_angle', 'angle_var', 'mean_to_angle', ...
         'mean_dest_angle', 'mean_from_angle', 'to_angle_var', 'dest_angle_var', ...
         'from_angle_var', 'mean_aux', 'aux_sd', 'mean_to_aux', 'mean_dest_aux', ...
         'mean_from_aux', 'to_aux_sd', 'dest_aux_sd', 'from_aux_sd'}));
end
