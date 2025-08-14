%% Mayar Ariss (ariss@mit.edu) - June 2025
%% Plot Corrected ADRS Curves for Different Thicknesses (Full Nodal Masses) + Drift Thresholds
clear; clc; close all;

modelFolder = 'C:/Users/mayar/OneDrive - Massachusetts Institute of Technology/FEM_buildings/results/Approach_B/p4_08_Country_house_3/EFM_analysis';
thicknesses = [0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50];
removeLastN = 10;

g = 9.81; % m/s²
H = 3.52;  % Wall height in meters

% Drift thresholds (for FEMA commentary)
driftLevels = [0.003, 0.006, 0.01]; % 0.3%, 0.6%, 1%
driftLabels = {'0.3\% Drift', '0.6\% Drift', '1\% Drift'};
dispThresholds = driftLevels * H;

% Nodal vertical masses per thickness [kg]
Mz_all = [
    1817, 431, 4585, 1210, 714, 1769;
    2044, 485, 5158, 1362, 803, 1991;
    2271, 539, 5731, 1513, 892, 2212;
    2498, 592, 6304, 1664, 982, 2433;
    2726, 646, 6877, 1815, 1071, 2654;
    2953, 700, 7450, 1967, 1160, 2875;
    3180, 754, 8023, 2118, 1249, 3097;
    3407, 808, 8596, 2269, 1338, 3318;
    3634, 862, 9169, 2421, 1428, 3539;
    3861, 916, 9743, 2572, 1517, 3760;
    4088, 969,10316, 2723, 1606, 3981;
    4316,1023,10889, 2874, 1695, 4202;
    4543,1077,11462, 3026, 1785, 4424;
];

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
    displ = dispData(:, 2);
    baseShear = sum(reacData(:, 2:end), 2);
    n = min(length(displ), length(baseShear));
    cutoff = max(1, n - removeLastN);
    displ = displ(1:cutoff);
    baseShear = baseShear(1:cutoff);

    M_eff = sum(Mz_all(i, :));
    Sa_g = (-baseShear) ./ (M_eff * g);

    plot(displ, Sa_g, '.-', 'LineWidth', 1.4, ...
        'Color', colors(mod(i-1, size(colors,1)) + 1, :), ...
        'DisplayName', sprintf('$t = %.3f\\, \\mathrm{m}$', t));

    % Interpolate and report Sa at drift thresholds
    for d = 1:length(dispThresholds)
        Sa_at_drift = interp1(displ, Sa_g, dispThresholds(d), 'linear', NaN);
        if ~isnan(Sa_at_drift)
            fprintf('Thickness %.3f m | Drift %.1f%% (%.1f mm) → Sa ≈ %.3f g\n', ...
                t, driftLevels(d)*100, dispThresholds(d)*1000, Sa_at_drift);
        end
    end
end

% Add vertical drift threshold lines
% for d = 1:length(dispThresholds)
%     xline(dispThresholds(d), '--k', driftLabels{d}, ...
%         'Interpreter', 'latex', 'LabelVerticalAlignment', 'bottom', 'FontSize', 9);
% end

xlabel('Control Node Displacement [$\mathrm{m}$]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Spectral Acceleration [$g$]', 'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');
grid on; box on;

% xlim([0, 0.05]);
% ylim([0, 2]);

outputFolder = fullfile("C:\Users\mayar\OneDrive - Massachusetts Institute of Technology\FEM_buildings\results\Approach_B\p4_08_Country_house_3\pushover_images", 'figures');
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

exportgraphics(gcf, fullfile(outputFolder, 'adrs_curves_thicknesses_drift.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(outputFolder, 'adrs_curves_thicknesses_drift.png'), 'Resolution', 600);
print(gcf, fullfile(outputFolder, 'adrs_curves_thicknesses_drift.svg'), '-dsvg', '-painters');

