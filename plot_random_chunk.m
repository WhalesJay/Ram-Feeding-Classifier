function plot_random_chunk(whale_names, chunk_length)
    % plot_random_chunk - Randomly selects and plots a chunk of data from a whale object.
    %
    % Parameters:
    % whale_names: List of whale object names.
    % chunk_length: Desired length of the continuous chunk where MouthOpen == 1.

    % Retry loop to find a valid whale object
    max_attempts = 10; % Limit the number of retries
    for attempt = 1:max_attempts
        random_index = randi(length(whale_names));
        random_whale_name = whale_names(random_index).name;

        % Check if the variable exists in the base workspace
        if evalin('base', sprintf('exist(''%s'', ''var'')', random_whale_name))
            % Evaluate the random whale object from the base workspace
            random_whale_obj = evalin('base', random_whale_name);
            break; % Exit the loop once a valid object is found
        end
    end

    % If no valid object is found after max_attempts
    if ~exist('random_whale_obj', 'var')
        error("No valid whale object found after %d attempts.", max_attempts);
    end

    % Extract necessary data
    datenum_values = random_whale_obj.Datenum;         % Time data
    %datenum_values = datetime(datenum_values, 'ConvertFrom', 'datenum');
    depth_values = random_whale_obj.Depth;            % Depth data
    mouth_open = random_whale_obj.MouthOpen;          % MouthOpen (1 or 0)
    speed = random_whale_obj.speed;                   % Speed data
    fluking_signal = random_whale_obj.fluking_signal; % Fluking signal data

    % Identify continuous chunks where MouthOpen == 1
    open_indices = find(mouth_open == 1); % Indices where MouthOpen == 1
    diff_indices = diff(open_indices);   % Difference between consecutive indices

    % Find continuous segments
    breaks = [0; find(diff_indices > 1); length(open_indices)];
    segments = arrayfun(@(x) open_indices(breaks(x)+1:breaks(x+1)), 1:length(breaks)-1, 'UniformOutput', false);

    % Filter segments with at least chunk_length rows
    valid_segments = segments(cellfun(@length, segments) >= chunk_length);

    % Ensure there are valid segments
    if isempty(valid_segments)
        error("No continuous chunks with at least %d rows where MouthOpen == 1.", chunk_length);
    end

    % Randomly select one valid segment
    selected_segment_open = valid_segments{randi(length(valid_segments))};

    % Select the first chunk_length rows from the chosen segment
    selected_indices_open = selected_segment_open(1:chunk_length);

    % Identify a random segment where MouthOpen == 0 with continuous indices
    closed_indices = find(mouth_open == 0);
    if isempty(closed_indices) || length(closed_indices) < chunk_length
        error("Not enough continuous rows with MouthOpen == 0 for the specified chunk length.");
    end

    % Find continuous segments for MouthClosed
    diff_closed_indices = diff(closed_indices);   % Difference between consecutive indices for MouthClosed
    closed_breaks = [0; find(diff_closed_indices > 1); length(closed_indices)];
    closed_segments = arrayfun(@(x) closed_indices(closed_breaks(x)+1:closed_breaks(x+1)), 1:length(closed_breaks)-1, 'UniformOutput', false);
    
    % Filter segments with at least chunk_length rows for MouthClosed
    valid_closed_segments = closed_segments(cellfun(@length, closed_segments) >= chunk_length);

    % Ensure there are valid continuous segments for MouthClosed
    if isempty(valid_closed_segments)
        error("No continuous chunks with at least %d rows where MouthOpen == 0.", chunk_length);
    end

    % Randomly select one valid segment for MouthClosed
    selected_segment_closed = valid_closed_segments{randi(length(valid_closed_segments))};

    % Select the first chunk_length rows from the chosen segment for MouthClosed
    selected_indices_closed = selected_segment_closed(1:chunk_length);

    % Extract the corresponding data for MouthOpen == 1
    datenum_open = datenum_values(selected_indices_open);
    depth_open = depth_values(selected_indices_open);
    speed_open = speed(selected_indices_open);
    fluking_open = fluking_signal(selected_indices_open);

    % Extract the corresponding data for MouthOpen == 0
    datenum_closed = datenum_values(selected_indices_closed);
    depth_closed = depth_values(selected_indices_closed);
    speed_closed = speed(selected_indices_closed);
    fluking_closed = fluking_signal(selected_indices_closed);

    % Plot the data
    figure;

    % Column 1: Mouth Open
    subplot(3, 2, 1); % Fluking Signal
    plot(datenum_open, fluking_open, '-o', 'Color', [0.8500, 0.3250, 0.0980]); % Orange
    xlabel('Time');
    ylabel('Fluking Signal');
    ylim([-1 1]);
    title('Mouth Open');

    subplot(3, 2, 3); % Speed
    plot(datenum_open, speed_open, '-o', 'Color', [0.4660, 0.6740, 0.1880]); % Green
    xlabel('Time');
    ylabel('Speed (m/s)');
    ylim([0 6]);

    subplot(3, 2, 5); % Depth
    plot(datenum_open, depth_open, '-o', 'Color', [0.4940, 0.1840, 0.5560]); % Purple
    xlabel('Time');
    ylabel('Depth (m)');
    ylim([0 80]);
    set(gca, 'YDir', 'reverse'); % Invert y-axis for depth

    % Column 2: Mouth Closed
    subplot(3, 2, 2); % Fluking Signal
    plot(datenum_closed, fluking_closed, '-o', 'Color',[0.8500, 0.3250, 0.0980]); % Orange
    xlabel('Time');
    ylabel('Fluking Signal');
    ylim([-1 1]);
    title('Mouth Closed');

    subplot(3, 2, 4); % Speed
    plot(datenum_closed, speed_closed, '-o', 'Color',[0.4660, 0.6740, 0.1880]); % Green
    xlabel('Time');
    ylabel('Speed (m/s)');
    ylim([0 6]);

    subplot(3, 2, 6); % Depth
    plot(datenum_closed, depth_closed, '-o', 'Color', [0.4940, 0.1840, 0.5560]); % Purple
    xlabel('Time');
    ylabel('Depth (m)');
    ylim([0 80]);
    set(gca, 'YDir', 'reverse'); % Invert y-axis for depth

    % Remove grid lines
    set(findall(gcf, 'Type', 'axes'), 'Layer', 'top');
    disp("Random chunks plotted for whale: " + random_whale_name);
end
