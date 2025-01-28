% Input coordinates
coordinates = {
    0 0 1 '40099';
    1 0 2 '400105';
    1 0 3 '400106';
    1 1 4 '500127';
    1 1 5 '500128';
    0 1 6 '500123';
    0 1 7 '500124';
    0 0 8 '400100';
    1 0 9 '10015';
    2 0 10 '10019';
    2 0 11 '10020';
    2 1 12 '20053';
    2 1 13 '20054';
    1 1 14 '20045';
    1 1 15 '20046';
    1 0 16 '10016'
};

% Extract the unique coordinates
[unique_coords, ~, idx] = unique(cell2mat(coordinates(:, 1:2)), 'rows', 'stable');

% Find groups of nodes sharing the same coordinates
groups = accumarray(idx, 1:numel(idx), [], @(x) {sort(x)});

% Iterate through groups
for i = 1:numel(groups)
    group_indices = groups{i};
    if numel(group_indices) == 4 % Check if group has four nodes
        % Append node names as specified
        coordinates{group_indices(3), 4} = [coordinates{group_indices(2), 4}];
    elseif numel(groups{i}) == 6 % Check if the group has six nodes
        % Replace node names for groups of six
        coordinates{groups{i}(2), 4} = coordinates{groups{i}(3), 4}; % Replace first with third
        coordinates{groups{i}(4), 4} = coordinates{groups{i}(5), 4}; % Replace second with fourth
    elseif numel(groups{i}) == 8 % Check if the group has eight nodes
        % Replace node names for groups of eight
        coordinates{groups{i}(1), 4} = coordinates{groups{i}(8), 4}; % Replace first with eighth
        coordinates{groups{i}(2), 4} = coordinates{groups{i}(3), 4}; % Replace second with third
        coordinates{groups{i}(4), 4} = coordinates{groups{i}(5), 4}; % Replace fourth with fifth
        coordinates{groups{i}(6), 4} = coordinates{groups{i}(7), 4}; % Replace sixth with seventh
    end
end

unique(coordinates(:,4))
% Display the modified coordinates
coordinates
