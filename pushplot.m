% Mayar Ariss (ariss@mit.edu) - June 2025
%% Plot Pushover Curves for Different Thicknesses (Filtered) + Drift Thresholds
clear; clc; close all;

modelFolder = 'C:/Users/mayar/OneDrive - Massachusetts Institute of Technology/FEM_buildings/results/Approach_B/p4_08_Country_house_3/EFM_analysis';
thicknesses = [0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50];
removeLastN = 10;

% Wall height
H = 3.52; % meters

% Drift thresholds non-infill (FEMA damage interpretation)
driftLevels = [0.003, 0.006, 0.01]; % 0.3%, 0.6%, 1%
driftLabels = {'0.1\% Drift', '0.3\% Drift', '0.6\% Drift'};
dispThresholds = driftLevels * H;

colors = [
    0.1216, 0.4667, 0.7059;
    1.0000, 0.4980, 0.0549;
    0.1725, 0.6275, 0.1725;
    0.8392, 0.1529, 0.1569;
    0.5804, 0.4039, 0.7412;
    0.5490, 0.3373, 0.2941;
    0.8902, 0.4667, 0.7608;
    0.4980, 0.4980, 0.4980;
    0.7373, 0.7412, 0.1333;
    0.0902, 0.7451, 0.8118;
    0.7373, 0.2353, 0.2392;
    0.3882, 0.6706, 0.8000;
    0.9686, 0.7137, 0.8235
];

figure('Color', 'w'); hold on;

for i = 1:length(thicknesses)
    t = thicknesses(i);
    suffix3 = sprintf('_%.3f', t);
    suffix2 = sprintf('_%.2f', t);

    dispFile3 = fullfile(modelFolder, 'outputFiles', ['output_disp', suffix3,'v2.out']);
    reacFile3 = fullfile(modelFolder, 'outputFiles', ['output_reac_all', suffix3,'v2.out']);
    dispFile2 = fullfile(modelFolder, 'outputFiles', ['output_disp', suffix2,'v2.out']);
    reacFile2 = fullfile(modelFolder, 'outputFiles', ['output_reac_all', suffix2,'v2.out']);

    if isfile(dispFile3) && isfile(reacFile3)
        dispFile = dispFile3; reacFile = reacFile3;
    elseif isfile(dispFile2) && isfile(reacFile2)
        dispFile = dispFile2; reacFile = reacFile2;
    else
        warning('Missing data for thickness %.3f m', t); continue;
    end

    dispData = load(dispFile);
    reacData = load(reacFile);

    % Check size
    if size(dispData, 2) < 2
        error('`%s` does not have at least 2 columns. Check file format or header.', dispFile);
    end

    displ = dispData(:, 2);
    baseShear = sum(reacData(:, 2:end), 2);

    n = min(length(displ), length(baseShear));
    cutoff = max(1, n - removeLastN);

    displ = displ(1:cutoff);
    baseShear = -baseShear(1:cutoff); % flip sign to match plotting convention

    % Plot pushover curve
    plot(displ, baseShear, '.-', 'LineWidth', 1.2, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('$t = %.3f\\, \\mathrm{m}$', t));

    % Interpolate and report base shear at drift thresholds
    for d = 1:length(dispThresholds)
        shear_at_drift = interp1(displ, baseShear, dispThresholds(d), 'linear', NaN);
        if ~isnan(shear_at_drift)
            fprintf('Thickness %.3f m | Drift %.1f%% (%.1f mm) → Base Shear ≈ %.1f N\n', ...
                t, driftLevels(d)*100, dispThresholds(d)*1000, shear_at_drift);
        end
    end
end

% Add vertical drift threshold lines
% for d = 1:length(dispThresholds)
%     xline(dispThresholds(d), '--k', driftLabels{d}, ...
%         'Interpreter', 'latex', 'LabelVerticalAlignment', 'bottom', 'FontSize', 9);
% end

xlabel('Control Node Displacement [$\mathrm{m}$]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Base Shear [$\mathrm{N}$]', 'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');
grid on; box on;

xlim([0, 0.01]);

% Save figure
outputFolder = fullfile("C:\Users\mayar\OneDrive - Massachusetts Institute of Technology\FEM_buildings\results\Approach_B\p4_08_Country_house_3\pushover_images", 'figures');
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

saveas(gcf, fullfile(outputFolder, 'pushover_curves_thicknesses_drift.pdf'));
saveas(gcf, fullfile(outputFolder, 'pushover_curves_thicknesses_drift.png'));
print(fullfile(outputFolder, 'pushover_curves_thicknesses_drift.svg'), '-dsvg', '-painters');
