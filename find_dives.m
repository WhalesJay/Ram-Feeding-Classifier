function dives = find_dives(p, mindepth, sampling_rate, surface, findall)
    % FIND_DIVES: Find the start and end time cues for dives or flights in a time series.
    % 
    % Inputs:
    %   - p: Depth or altitude time series in meters (vector or struct with fields `data` and `sampling_rate`).
    %   - mindepth: Threshold in meters to recognize a dive or flight.
    %   - sampling_rate: Sampling rate of the time series in Hz.
    %   - surface: Threshold in meters where the animal is presumed to be at the surface (default = 1).
    %   - findall: Boolean to include incomplete dives at start/end of the record (default = false).
    % 
    % Output:
    %   - dives: A table with columns:
    %       - start: Start time (seconds) of each dive/flight.
    %       - end: End time (seconds) of each dive/flight.
    %       - max: Maximum depth/altitude reached during the dive/flight.
    %       - tmax: Time (seconds) at which max depth/altitude occurred.

    if nargin < 2
        error("Inputs for 'p' and 'mindepth' are required.");
    end
    if isstruct(p)
        if isfield(p, 'data') && isfield(p, 'sampling_rate')
            sampling_rate = p.sampling_rate;
            p = p.data;
        else
            error("If 'p' is a struct, it must contain 'data' and 'sampling_rate' fields.");
        end
    end

    if nargin < 4 || isempty(surface)
        surface = 1; % Default surface threshold
    end
    if nargin < 5 || isempty(findall)
        findall = false; % Default findall flag
    end

    searchlen = 20; % Lookback/forward time in seconds to find surface
    dpthresh = 0.25; % Threshold for vertical velocity to detect surface
    dp_lp = 0.25; % Low-pass filter frequency for vertical velocity

    % Ensure p is a column vector
    if isrow(p)
        p = p';
    end

    % Adjust first depth value if it starts deeper than mindepth
    if p(1) > mindepth
        p(1) = mindepth - 0.25;
    end

    % Detect threshold crossings and surface times
    tth = find(diff(p > mindepth) > 0); % Start of dives
    tsurf = find(p < surface); % Surface points
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
                    toff(k) = length(p);
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
    dp = filter(fir1(n, dp_lp / (sampling_rate / 2)), 1, [0; diff(p)] * sampling_rate);

    dmax = zeros(length(ton), 2);
    for k = 1:length(ton)
        ind_start = max(ton(k) - round(searchlen * sampling_rate), 1):ton(k);
        ind_end = toff(k):min(toff(k) + round(searchlen * sampling_rate), length(p));

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

        [dm, km] = max(p(ton(k):toff(k)));
        dmax(k, :) = [dm, (ton(k) + km - 1) / sampling_rate];
    end

    % Compile results into a table
    start_times = ton / sampling_rate;
    end_times = toff / sampling_rate;
    dives = table(start_times, end_times, dmax(:, 1), dmax(:, 2), ...
                  'VariableNames', {'start', 'end', 'max', 'tmax'});
end
