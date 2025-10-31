% Extract data
%Change this to select whatever you want to select
%figure_data = rw24070642processed(220000:230000, :); 

%Classified deep dive
%figure_data = rw24070869processed(29000:35000, :);

%Classified Shallow Dive
figure_data = rw24070642processed(220000:230000, :); 

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
%ylabel('Fluking Signal (rad)');
set(gca, 'XTickLabel', []); % Remove x-axis labels for clarity
set(gca, 'XTick', []); % Remove x-axis ticks
set(gca, 'FontSize', 20); % Adjust 14 to desired size

hold off;

%% Plot Speed (Middle, 1/4 height)
subplot(5,1,2);
hold on;
colors = {'#B1CEDF', '#06559F'}; % {Colombia Blue, Lapis Lazuli}
for i = 1:length(fig_time) - 1
    line_color = colors{fig_open(i) + 1}; % Select color based on MouthOpen status
    plot(fig_time(i:i+1), fig_speed(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
%ylabel('Speed (m/s)');
set(gca, 'XTickLabel', []); % Remove x-axis labels for clarity
set(gca, 'XTick', []); % Remove x-axis ticks
set(gca, 'FontSize', 20); % Adjust 14 to desired size

hold off;

%% Plot Depth (Bottom, 1/2 height)
subplot(5,1,[3 5]); % Make depth take half of the total figure height
hold on;
colors = {'#D2E4C4', '#4E8056'}; % {Tea Green, Fern Green}
for i = 1:length(fig_time) - 1
    line_color = colors{fig_open(i) + 1}; % Select color based on MouthOpen status
    plot(fig_time(i:i+1), fig_depth(i:i+1), '-', 'Color', line_color, 'LineWidth', 1.5);
end
ylim([0 50]);
set(gca, 'YDir', 'reverse'); % Reverse y-axis for depth
%ylabel('Depth (m)');
xlabel('Time (ADT)');
set(gca, 'FontSize', 20); % Adjust 14 to desired size

hold off;

