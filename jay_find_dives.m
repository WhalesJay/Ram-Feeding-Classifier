function data = jay_find_dives(data, mindepth, sampling_rate, surface, findall, na_rm, prop)
    % JAY_FIND_DIVES: Combines dive detection, statistical analysis, and classification into one function.
    %
    % Inputs:
    %   - data: Struct containing trial data with fields `Depth` (depth data) and `pitch` (angular data).
    %   - mindepth: Threshold in meters to recognize a dive or flight.
    %   - sampling_rate: Sampling rate of the time series in Hz.
    %   - surface: Threshold in meters where the animal is presumed to be at the surface (default = 1).
    %   - findall: Boolean to include incomplete dives at start/end of the record (default = false).
    %   - na_rm: Boolean to indicate whether to remove NaN values from auxiliary data (default = true).
    %   - prop: proportion of time spent at bottom, deafault 0.85

    % Outputs:
    %   - Dive statistics and classification are added to the data struct, and a plot is generated.

    % Step 1: Detect dives
    if nargin < 3
        error("Not enough input arguments for dive detection.");
    end

    if nargin < 4 || isempty(surface)
        surface = 1; % Default surface threshold
    end
    if nargin < 5 || isempty(findall)
        findall = false; % Default findall flag
    end
    
    if nargin < 6 || isempty(na_rm)
        na_rm = true;
    end
    if nargin < 7 || isempty(prop)
        prop = 0.85;
    end
    searchlen = 20; % Lookback/forward time in seconds to find surface
    dpthresh = 0.25; % Threshold for vertical velocity to detect surface
    dp_lp = 0.25; % Low-pass filter frequency for vertical velocity

    % Ensure data.Depth is a column vector
    if isrow(data.Depth)
        data.Depth = data.Depth';
    end

    % Adjust first depth value if it starts deeper than mindepth
    if data.Depth(1) > mindepth
        data.Depth(1) = mindepth - 0.25;
    end

    % Detect threshold crossings and surface times
    tth = find(diff(data.Depth > mindepth) > 0); % Start of dives
    tsurf = find(data.Depth < surface); % Surface points
    ton = zeros(size(tth)); % Dive start times
    toff = ton; % Dive end times

    k = 0;
    for kth = 1:length(tth)
        if all(tth(kth) > toff)
            ks0 = tsurf(tsurf < tth(kth));
            ks1 = tsurf(tsurf > tth(kth));

            if findall || (~isempty(ks0) && ~isempty(ks1))
                k = k + 1;

                if isempty(ks0)
                    ton(k) = 1;
                else
                    ton(k) = max(ks0);
                end

                if isempty(ks1)
                    toff(k) = length(data.Depth);
                else
                    toff(k) = min(ks1);
                end
            end
        end
    end

    % Truncate ton and toff arrays
    ton = ton(1:k);
    toff = toff(1:k);

    % Compute vertical velocity and refine surface times
    n = round(4 * sampling_rate / dp_lp);
    dp = filter(fir1(n, dp_lp / (sampling_rate / 2)), 1, [0; diff(data.Depth)] * sampling_rate);

    dmax = zeros(length(ton), 2);
    for k = 1:length(ton)
        ind_start = max(ton(k) - round(searchlen * sampling_rate), 1):ton(k);
        ind_end = toff(k):min(toff(k) + round(searchlen * sampling_rate), length(data.Depth));

        ki_start = find(dp(ind_start) < dpthresh, 1, 'last');
        ki_end = find(dp(ind_end) > -dpthresh, 1, 'first');

        if isempty(ki_start)
            ki_start = 1; % Default to the first index if not found
        end
        if isempty(ki_end)
            ki_end = length(ind_end); % Default to the last index if not found
        end

        ton(k) = ind_start(max(1, ki_start));
        toff(k) = ind_end(min(length(ind_end), ki_end));

        [dm, km] = max(data.Depth(ton(k):toff(k)));
        dmax(k, :) = [dm, (ton(k) + km - 1) / sampling_rate];
    end

    % Create the dive cues table
    start_times = ton / sampling_rate;
    end_times = toff / sampling_rate;
    dive_cues = table(start_times, end_times, dmax(:, 1), dmax(:, 2), ...
                      'VariableNames', {'start', 'end', 'max', 'tmax'});

    % Step 2: Calculate dive statistics

    % Convert dive_cues from table to numeric array and multiply by sampling rate
    start_times = table2array(dive_cues(:, 1));  % Extract start times from the first column
    end_times = table2array(dive_cues(:, 2));    % Extract end times from the second column
    di = round([start_times, end_times] * sampling_rate);   % Multiply by sampling rate and round the values

    % Initialize output dataframe
    Y = table((1:size(dive_cues, 1))', 'VariableNames', {'num'});

    % Loop through each dive/flight
    for d = 1:size(dive_cues, 1)
        z = data.Depth(di(d, 1):di(d, 2));  % Get depth data for the current dive
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
        
        % Process angular data (pitch)
        a = data.pitch(di(d, 1):di(d, 2));  % Get angular data for the current dive
        at = a(1:pt(1));  % Data for the "to" phase
        af = a(pt(end):end);  % Data for the "from" phase
        ad = a(pt(1):pt(end));  % Data for the "destination" phase
        Y.st(d) = start_times(d); 
        Y.et(d) = end_times(d);
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
    end
    for i = 1:height(Y)
    Y.bottom_dur_perc(i) = Y.dest_dur(i)/Y.dur(i);
    end

    Y.diveshape = repmat("unknown", height(Y), 1); % Initialize with default value


    Y.diveshape(Y.bottom_dur_perc > 0.5) = "square";
    Y.diveshape(Y.bottom_dur_perc >= 0.2 & Y.bottom_dur_perc <= 0.5) = "U";
    Y.diveshape(Y.bottom_dur_perc < 0.2) = "V";

    % Step 3: Classify and plot the data
    data.class = repmat({'sur'}, height(data), 1); % Default classification is 'sur'
    data.shape = repmat("unknown", height(data), 1); % Initialize with default value
    data.reltime = (1:height(data))' / 10;
    % Loop through the dive statistics table Y to classify data rows
    for d = 1:height(Y)
        % Find the rows in data that correspond to the current dive/flight
        trial_idx = data.reltime >= Y.st(d) & ...
                     data.reltime <= Y.et(d);
        data.shape(trial_idx) = Y.diveshape(d);

        
        % Mark as 'dest' for destination phase
        dest_idx = data.reltime >= Y.dest_st(d) & ...
                   data.reltime <= Y.dest_et(d);
        data.class(dest_idx) = {'dest'};
        
        % Mark as 'to' for "to" phase
        to_idx = data.reltime >= Y.st(d) & ...
                 data.reltime <= Y.dest_st(d);
        data.class(to_idx) = {'to'};
        
        % Mark as 'from' for "from" phase
        from_idx = data.reltime >= Y.dest_et(d)& ...
                   data.reltime <= Y.et(d);
        data.class(from_idx) = {'from'};
    end

    % Plot the classified data with inverted y-axis
    figure;
    gscatter(data.Datenum, data.Depth, data.class);
    ylabel('Depth');
    xlabel('Time');
    title('Dive Phase Classification');
    set(gca, 'YDir', 'reverse'); % Invert the Y-axis
    legend('show');
end
