% Mayar Ariss (ariss@mit.edu) - June 2025
%% Setup
clear; clc; close all;
set(groot,'defaultAxesTickLabelInterpreter','default');
modelPath = 'C:/Users/mayar/OneDrive - Massachusetts Institute of Technology/FEM_buildings/results/Approach_B/p4_09_Country_house_4/EFM_analysis/results.mat';
outputFolder = fullfile("C:\Users\mayar\OneDrive - Massachusetts Institute of Technology\FEM_buildings\results\Approach_B\p4_08_Country_house_3\pushover_images", 'figures'); % <- create a 'figures' folder if needed

load(modelPath);

periods = [
     0.1738
    0.1605
    0.1125
    0.1114
    0.0878
    0.0831
    0.0713
    0.0658
    0.0529
    0.0476
];

nModes = 10;
nDOF_per_node = 6;
modeFolder = 'C:/Users/mayar/OneDrive - Massachusetts Institute of Technology/FEM_buildings/results/Approach_B/p4_08_Country_house_3/EFM_analysis/outputFiles/Modal1/';

wallsToAnalyze = [2,3];

% Storage
rms_OOP = zeros(nModes, length(wallsToAnalyze));
rms_IP = zeros(nModes, length(wallsToAnalyze));
modeShapesOOP = cell(nModes, length(wallsToAnalyze));
modeShapesIP = cell(nModes, length(wallsToAnalyze));

% --- Loop over walls and modes
for w = 1:length(wallsToAnalyze)
    targetWall = wallsToAnalyze(w);
    wall_idx = find(arrayfun(@(n) ~isempty(n.wall) && any(n.wall == targetWall), model.node));
    
    for modeNum = 1:nModes
        modeFile = sprintf('%s/EFM_mode%d.out', modeFolder, modeNum);
        modeData = load(modeFile);
        
        if mod(numel(modeData), nDOF_per_node) ~= 0
            error('Mode %d file corrupted.', modeNum);
        end
        nNodes = numel(modeData) / nDOF_per_node;
        modeMatrix = reshape(modeData, [nDOF_per_node, nNodes])';
        
        % Extract displacements
        UX = modeMatrix(wall_idx, 1);
        UY = modeMatrix(wall_idx, 2);
        UZ = modeMatrix(wall_idx, 3);
        
        % Store RMS for each type
        rms_OOP(modeNum, w) = sqrt(mean(UY.^2));                     
        rms_IP(modeNum, w)  = sqrt(mean((UX.^2 + UZ.^2)/2));          
        
        % Save mode shapes
        modeShapesOOP{modeNum, w} = UY;
        modeShapesIP{modeNum, w} = sqrt((UX.^2 + UZ.^2)/2);
    end
end


%% Figure 1: Mode Shapes for Wall 2 and Wall 3 (Separate Subplots, one global legend)
figure('Color','w');

selectedModes = 1:10;
% --- Define Custom Colors
colors = [
    0,    114, 178;   % Strong Blue
    213,  94,  0;     % Vermillion (Strong Red-Orange)
    0,    158, 115;   % Bluish Green
    204,  121, 167;   % Reddish Purple
    230,  159, 0;     % Orange
    86,   180, 233;   % Sky Blue
    240,  228, 66;    % Yellow
    0,    0,   204;   % Bluish Purple
    0,    128, 0;     % Olive Green
    64,   64,  64;    % Blackish Gray
] / 255;


% --- Subplot Wall 2
subplot(2,1,1); hold on;
for i = 1:length(selectedModes)
    modeNum = selectedModes(i);
    outOfPlane_disp = modeShapesOOP{modeNum, 1};
    p2(i) = plot(1:length(outOfPlane_disp), outOfPlane_disp*1000, '.-', 'LineWidth', 0.7, ...
        'Color', colors(i,:));
end
xlabel('Node Index', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('OOP Displacement [mm]', 'Interpreter', 'latex', 'FontSize', 10);
xlim([1, length(modeShapesOOP{modeNum, 1})])
title('Wall 2', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);
grid on;

% --- Subplot Wall 3
subplot(2,1,2); hold on;
for i = 1:length(selectedModes)
    modeNum = selectedModes(i);
    outOfPlane_disp = modeShapesOOP{modeNum, 2};
    plot(1:length(outOfPlane_disp), outOfPlane_disp*1000, '.-', 'LineWidth', 0.7, ...
        'Color', colors(i,:));
end
xlabel('Node Index', 'Interpreter', 'latex', 'FontSize', 10);
xlim([1, length(modeShapesOOP{modeNum, 2})])
ylabel('OOP Displacement [mm]', 'Interpreter', 'latex', 'FontSize', 10);
title('Wall 3', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);
grid on;

% --- Global legend
%sgtitle('Out-of-Plane Mode Shapes (Wall 2 and Wall 3)', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 16);
allAxes = findall(gcf, 'type', 'axes', '-not', 'Tag', 'legend');
for k = 1:length(allAxes)
    ax = allAxes(k);
    pos = get(ax, 'Position');
    pos(3) = pos(3) * 0.8;  % shrink width to 80%
    set(ax, 'Position', pos);
end
hBigAx = axes(gcf, 'Visible', 'off', 'Position', [0 0 1 1]);
legend(hBigAx, p2, arrayfun(@(x) sprintf('Mode %d (T=%.2fs)', x, periods(x)), selectedModes, 'UniformOutput', false), ...
    'Orientation', 'vertical', ...
    'Interpreter', 'latex', ...
    'Box', 'on', ...
    'FontSize', 7, ...
    'Location', 'eastoutside');

saveas(gcf, fullfile(outputFolder, 'OutOfPlane_ModeShapes.pdf'));
print(gcf, fullfile(outputFolder, 'OutOfPlane_ModeShapes.svg'), '-dsvg', '-painters');

%% Print numerical summaries
disp('RMS Out-of-Plane [m]:');
disp(rms_OOP);
disp('RMS In-Plane [m]:');
disp(rms_IP);

for w = 1:length(wallsToAnalyze)
    [max_rms_OOP, idx_OOP] = max(rms_OOP(:,w));
    [max_rms_IP, idx_IP] = max(rms_IP(:,w));
    fprintf('Wall %d: Max OOP RMS = %.4e m at Mode %d\n', wallsToAnalyze(w), max_rms_OOP, idx_OOP);
    fprintf('Wall %d: Max IP RMS = %.4e m at Mode %d\n', wallsToAnalyze(w), max_rms_IP, idx_IP);
end
disp('Periods [s] of the analyzed modes:');
disp(periods);

wallsToAnalyze = [1, 2, 3, 4];

% Storage
rms_OOP = zeros(nModes, length(wallsToAnalyze));
rms_IP  = zeros(nModes, length(wallsToAnalyze));

% --- Loop over walls and modes
for w = 1:length(wallsToAnalyze)
    targetWall = wallsToAnalyze(w);
    wall_idx = find(arrayfun(@(n) ~isempty(n.wall) && any(n.wall == targetWall), model.node));
    
    for modeNum = 1:nModes
        modeFile = sprintf('%s/EFM_mode%d.out', modeFolder, modeNum);
        modeData = load(modeFile);
        
        if mod(numel(modeData), nDOF_per_node) ~= 0
            error('Mode %d file corrupted.', modeNum);
        end
        nNodes = numel(modeData) / nDOF_per_node;
        modeMatrix = reshape(modeData, [nDOF_per_node, nNodes])';
        
        UX = modeMatrix(wall_idx, 1);
        UY = modeMatrix(wall_idx, 2);
        UZ = modeMatrix(wall_idx, 3);

        % Component selection
        if ismember(targetWall, [1, 3])
            OOP_disp = UX;
            IP_disp  = sqrt((UY.^2 + UZ.^2)/2);
        else
            OOP_disp = UY;
            IP_disp  = sqrt((UX.^2 + UZ.^2)/2);
        end

        rms_OOP(modeNum, w) = sqrt(mean(OOP_disp.^2));
        rms_IP(modeNum, w)  = sqrt(mean(IP_disp.^2));
    end
end

%% Figure 2: RMS vs Mode and RMS vs Period (All 4 Walls, one global legend)
figure('Color','w');

% --- Colors (same as before)
colors = [
    44, 160, 44;    % Wall 1: Medium Green
    31, 119, 180;   % Wall 2: Strong Blue
    255, 127, 14;   % Wall 3: Strong Orange
    214, 39, 40;    % Wall 4: Bold Red
] / 255;

% Upstream fix: reorder walls to plot as 4,2,3,1 but label as 1,2,3,4
wallsToAnalyze = [2,3];  % you want walls 2 and 3

% Select the correct columns:
wallsIdx = [2, 3];  % explicitly: column 2 and 3
rms_OOP_sel = rms_OOP(:, wallsIdx);
rms_IP_sel  = rms_IP(:,  wallsIdx);

legendEntries = {};
legendHandles = [];

% --- Subplot 1: RMS vs Mode
subplot(2,1,1); hold on;
for w = 1:length(wallsIdx)
    colIdx = wallsIdx(w);
    wallNum = wallsToAnalyze(w);
    %displayedWall = wallsToAnalyze(w);

    labelWall = wallNum;
    if wallNum  == 1
        labelWall = 4;
    elseif wallNum  == 4
        labelWall = 1;
    end

    % OOP
    h1 = plot(1:nModes, rms_OOP(:,colIdx)*1000, '-o', 'Color', colors(colIdx,:), 'LineWidth', 0.7, 'MarkerSize', 3);
    legendEntries{end+1} = sprintf('Wall %d - OOP', labelWall);
    legendHandles(end+1) = h1;
    % IP
    h2 = plot(1:nModes, rms_IP(:,colIdx)*1000, '--s', 'Color', colors(colIdx,:), 'LineWidth', 0.7, 'MarkerSize', 3);
    legendEntries{end+1} = sprintf('Wall %d - IP', labelWall);
    legendHandles(end+1) = h2;
end
xlabel('Mode Number', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('RMS Displacement [mm]', 'Interpreter', 'latex', 'FontSize', 10);
%title('RMS Displacement vs Mode (All Walls)', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);
grid on;
xlim([1 nModes]);

% --- Subplot 2: RMS vs Period
subplot(2,1,2); hold on;
for w = 1:length(wallsIdx)
    colIdx = wallsIdx(w);
    wallNum = wallsToAnalyze(w);
    
    labelWall = wallNum;
    if wallNum == 1
        labelWall = 4;
    elseif wallNum == 4
        labelWall = 1;
    end

    plot(periods, rms_OOP(:,colIdx)*1000, '-o', 'Color', colors(colIdx,:), 'LineWidth', 0.7, 'MarkerSize', 3);
    plot(periods, rms_IP(:,colIdx)*1000, '--s', 'Color', colors(colIdx,:), 'LineWidth', 0.7, 'MarkerSize', 3);
end
xlabel('Period [s]', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('RMS Displacement [mm]', 'Interpreter', 'latex', 'FontSize', 10);
%title('RMS Displacement vs Period (All Walls)', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);
grid on;
xlim([min(periods) max(periods)]);

% --- Adjust and add legend
allAxes = findall(gcf, 'type', 'axes', '-not', 'Tag', 'legend');
for k = 1:length(allAxes)
    pos = get(allAxes(k), 'Position');
    pos(3) = pos(3) * 0.8;
    set(allAxes(k), 'Position', pos);
end
hBigAx = axes(gcf, 'Visible', 'off', 'Position', [0 0 1 1]);
legend(hBigAx, legendHandles, legendEntries, ...
    'Orientation', 'vertical', ...
    'Interpreter', 'latex', ...
    'Box', 'on', ...
    'FontSize', 7, ...
    'Location', 'eastoutside');

% --- Save figure
saveas(gcf, fullfile(outputFolder, 'RMS_OOP_IP_Walls1to4_vs_Mode_Period.pdf'));
print(gcf, fullfile(outputFolder, 'RMS_OOP_IP_Walls1to4_vs_Mode_Period.svg'), '-dsvg', '-painters');


%% Print Numerical Summary
disp('Out-of-Plane RMS [m]:');
disp(rms_OOP);

disp('In-Plane RMS [m]:');
disp(rms_IP);

for w = 1:length(wallsToAnalyze)
    [max_rms_OOP, idx_OOP] = max(rms_OOP(:,w));
    [max_rms_IP,  idx_IP ] = max(rms_IP(:,w));
    fprintf('Wall %d: Max OOP RMS = %.4e m at Mode %d\n', wallsToAnalyze(w), max_rms_OOP, idx_OOP);
    fprintf('Wall %d: Max IP  RMS = %.4e m at Mode %d\n', wallsToAnalyze(w), max_rms_IP,  idx_IP );
end
disp('==== SUMMARY OF RESULTS PER WALL ====');
for w = 1:length(wallsToAnalyze)
    wallNum = wallsToAnalyze(w);
    [max_rms_OOP, idx_OOP] = max(rms_OOP(:,w));
    [max_rms_IP, idx_IP ] = max(rms_IP(:,w));
    fprintf('Wall %d:\n', wallNum);
    fprintf('  Max OOP RMS = %.4f mm at Mode %d\n', max_rms_OOP*1000, idx_OOP);
    fprintf('  Max IP  RMS = %.4f mm at Mode %d\n', max_rms_IP*1000, idx_IP);
    fprintf('  OOP vs IP ratio (max values) = %.2f\n', max_rms_OOP / max_rms_IP);
    fprintf('----------------------------------------\n');
end

fprintf('\nPeriods [s]:\n');
disp(periods');
[max_rms_OOP_wall3, idx_OOP_wall3] = max(rms_OOP(:, 3));
fprintf('Wall 3: Max OOP RMS = %.4f mm at Mode %d\n', max_rms_OOP_wall3 * 1000, idx_OOP_wall3);
