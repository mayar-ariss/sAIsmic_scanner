% Define your matrices XNodes_filtered, YNodes_filtered, XSpandrels_filtered,
% YSpandrels_filtered, XPiers_filtered, and YPiers_filtered

% Define a function to find values occurring four times consecutively
find_consecutive_four = @(matrix) unique(matrix(all(diff(matrix) == 0, 1), :), 'rows');

% Find X_face for each matrix
X_face = [];
X_matrices = {round(XNodes_filtered, 4), round(XSpandrels_filtered, 4), round(XPiers_filtered, 4)};
for i = 1:numel(X_matrices)
    matrix = X_matrices{i};
    [rows, cols] = size(matrix);
    for col = 1:cols
        for row = 1:rows-3
            if all(matrix(row:row+3, col) == matrix(row, col))
                X_face = [X_face; matrix(row, col)];
            end
        end
    end
end

% Extract unique values
X_face = unique(X_face);

% Find X_face for each matrix
Y_face = [];
Y_matrices = {round(YNodes_filtered, 4), round(YSpandrels_filtered, 4), round(YPiers_filtered, 4)};
for i = 1:numel(Y_matrices)
    matrix = Y_matrices{i};
    [rows, cols] = size(matrix);
    for col = 1:cols
        for row = 1:rows-3
            if all(matrix(row:row+3, col) == matrix(row, col))
                Y_face = [Y_face; matrix(row, col)];
            end
        end
    end
end

% Extract unique values
Y_face = unique(Y_face);


% Define your matrices XNodes_filtered, YNodes_filtered, XSpandrels_filtered,
% YSpandrels_filtered, XPiers_filtered, and YPiers_filtered

% Define a function to find values occurring four times consecutively
find_consecutive_four = @(matrix) unique(matrix(all(diff(matrix) == 0, 1), :));

% Initialize variables
X_face = [];
Y_face = [];
X_matrices = {round(XNodes_filtered, 4), round(XSpandrels_filtered, 4), round(XPiers_filtered, 4)};
Y_matrices = {round(YNodes_filtered, 4), round(YSpandrels_filtered, 4), round(YPiers_filtered, 4)};

% Find X_face for each matrix
for i = 1:numel(X_matrices)
    matrix = X_matrices{i};
    [rows, cols] = size(matrix);
    for col = 1:cols
        for row = 1:rows-3
            if all(matrix(row:row+3, col) == matrix(row, col))
                X_face = [X_face; matrix(row, col)];
            end
        end
    end
end

% Find Y_face for each matrix
for i = 1:numel(Y_matrices)
    matrix = Y_matrices{i};
    [rows, cols] = size(matrix);
    for col = 1:cols
        for row = 1:rows-3
            if all(matrix(row:row+3, col) == matrix(row, col))
                Y_face = [Y_face; matrix(row, col)];
            end
        end
    end
end

% Extract unique values
X_face = unique(X_face);
Y_face = unique(Y_face);


