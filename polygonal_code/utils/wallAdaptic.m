% Mayar Ariss 25 Mar 2024

function [successOutput, node_name] = wallAdaptic(XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face)

%[successOutput, node_name] = wallAdaptic(XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face)

%[successOutput, XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face]=generateModel(readTremuriInput('house1.txt'), 'ColorPiers', [0.2 0.2 1.0], 'ColorSpandrels', [0.2 1. .2], 'styleNodes', 'sk')

%FILE NAME
fn='_BUILDING2_nd';

% Combine X_face and Y_face into a single cell array to generate all 4
% walls
faces = {X_face, Y_face};

for face_idx = 1:length(faces)
    face_values = faces{face_idx};
    for i = 1:length(face_values)
        values = face_values(i);
        zero_indices = [];

        %% NODES
        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XNodes_filtered, 4) == values);
            prefix = 'W' + string(i) + '_X'+fn; %for when files are saved below
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YNodes_filtered, 4) == values);
            prefix = 'W' + string(i) + '_Y'+fn; %for when files are saved below
        end

        % Select X, Y, and Z for only front frame
        X_nodes = round(XNodes_filtered(:, zero_indices), 3);
        Y_nodes = round(YNodes_filtered(:, zero_indices), 3);
        Z_nodes = round(ZNodes_filtered(:, zero_indices), 3);

        % Convert to vectors
        X_nodes = X_nodes(:);
        Y_nodes = Y_nodes(:);
        Z_nodes = Z_nodes(:);

        % Rearrange every four-element cluster
        for i = 1:8:numel(X_nodes)
            % Extract the current four-element cluster
            cluster_x = X_nodes(i:i+3);
            cluster_y = Y_nodes(i:i+3);
            cluster_z = Z_nodes(i:i+3);

            % Sort the cluster
            sorted_cluster_x = sort(cluster_x);
            sorted_cluster_y = sort(cluster_y);
            sorted_cluster_z = sort(cluster_z);

            % Rearrange the cluster as a,b,b,a where a<b
            X_nodes(i:i+3) = sorted_cluster_x([2, 3, 3, 3]);
            X_nodes(i+4:i+7) = sorted_cluster_x([3, 2, 2, 2]);

            Y_nodes(i:i+3) = sorted_cluster_y([2, 3, 3, 3]);
            Y_nodes(i+4:i+7) = sorted_cluster_y([3, 2, 2, 2]);

            Z_nodes(i:i+3) = sorted_cluster_z([2, 2, 2, 3]);
            Z_nodes(i+4:i+7) = sorted_cluster_z([3, 3, 3, 2]);

        end


        %% SPANDRELS

        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XSpandrels_filtered, 4) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YSpandrels_filtered, 4) == values);
        end

        % Select X, Y, and Z for only front frame
        X_spandrels = round(XSpandrels_filtered(:, zero_indices), 3);
        Y_spandrels = round(YSpandrels_filtered(:, zero_indices), 3);
        Z_spandrels = round(ZSpandrels_filtered(:, zero_indices), 3);

        % Reshape X, Y, and Z matrices to column vectors
        X_spandrels = X_spandrels(:);
        Y_spandrels = Y_spandrels(:);
        Z_spandrels = Z_spandrels(:);

        % Rearrange every four-element cluster
        for i = 1:8:numel(X_spandrels)
            % Extract the current four-element cluster
            cluster_x = X_spandrels(i:i+3);
            cluster_y = Y_spandrels(i:i+3);
            cluster_z = Z_spandrels(i:i+3);

            % Sort the cluster
            sorted_cluster_x = sort(cluster_x);
            sorted_cluster_y = sort(cluster_y);
            sorted_cluster_z = sort(cluster_z);

            % Rearrange the cluster as a,b,b,a where a<b
            X_spandrels(i:i+3) = sorted_cluster_x([2, 3, 3, 3]);
            X_spandrels(i+4:i+7) = sorted_cluster_x([3, 2, 2, 2]);

            Y_spandrels(i:i+3) = sorted_cluster_y([2, 3, 3, 3]);
            Y_spandrels(i+4:i+7) = sorted_cluster_y([3, 2, 2, 2]);

            Z_spandrels(i:i+3) = sorted_cluster_z([2, 2, 2, 3]);
            Z_spandrels(i+4:i+7) = sorted_cluster_z([3, 3, 3, 2]);
        end


        %% PIERS

        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XPiers_filtered, 4) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YPiers_filtered, 4) == values);
        end


        X_piers = round(XPiers_filtered(:, zero_indices), 3);
        Y_piers = round(YPiers_filtered(:, zero_indices), 3);
        Z_piers = round(ZPiers_filtered(:, zero_indices), 3);

        % Reshape XPiers, YPiers, and ZPiers matrices to column vectors
        X_piers = X_piers(:);
        Y_piers = Y_piers(:);
        Z_piers = Z_piers(:);

        % Rearrange every four-element cluster
        for i = 1:8:numel(X_piers)
            % Extract the current four-element cluster
            cluster_x = X_piers(i:i+3);
            cluster_y = Y_piers(i:i+3);
            cluster_z = Z_piers(i:i+3);

            % Sort the cluster
            sorted_cluster_x = sort(cluster_x);
            sorted_cluster_y = sort(cluster_y);
            sorted_cluster_z = sort(cluster_z);

            % Rearrange the cluster as a,b,b,a where a<b
            X_piers(i:i+3) = sorted_cluster_x([2, 3, 3, 3]);
            X_piers(i+4:i+7) = sorted_cluster_x([3, 2, 2, 2]);

            Y_piers(i:i+3) = sorted_cluster_y([2, 3, 3, 3]);
            Y_piers(i+4:i+7) = sorted_cluster_y([3, 2, 2, 2]);

            Z_piers(i:i+3) = sorted_cluster_z([2, 2, 2, 3]);
            Z_piers(i+4:i+7) = sorted_cluster_z([3, 3, 3, 2]);
        end


        %% ADD RESTRAINTS

        % Determine the number of nodes, spandrels, piers and elements
        num_nodes= numel(X_nodes);
        num_spandrels = numel(X_spandrels);
        num_piers = numel(X_piers);
        num_elements=num_spandrels+num_piers+num_nodes;

        % Create a column vector to indicate whether each element is a spandrel (1) or not (0)
        is_spandrel = zeros(num_spandrels, 1); % Assuming all elements are spandrels initially
        % Create a combined matrix of X, Y, Z, their indices, and the spandrel indicator
        spandrels = [X_spandrels, Y_spandrels, Z_spandrels, (1:num_spandrels)', is_spandrel];

        patterns = {'ry+rz', 'rx+ry'};
        % Convert data to a cell array
        spandrels = num2cell(spandrels);
        % Initialize the counter for the pattern index
        pattern_index = 1;
        % Loop through the rows of data
        for i = 1:size(spandrels, 1)
            % Get the pattern for the current row
            pattern = patterns{pattern_index};
            % Assign the pattern to the 6th column of the current row
            spandrels{i, 6} = pattern;
            if mod(i, 4) == 0
                % Update the pattern index for every 4 rows
                pattern_index = mod(pattern_index, numel(patterns)) + 1;
            end
        end
        % Sort the spandrel combined matrix based on Z first, then Y, and finally on X within each Z-Y group
        sorted_spandrels = sortrows(spandrels, [3, 2, 1]);

        % Create a column vector to indicate whether each element is a piers (1) or not (0)
        is_pier = ones(num_piers, 1); % Assuming all elements are not piers initially
        % Create a combined matrix of XPiers, YPiers, ZPiers, their indices, and the piers indicator
        piers = [X_piers, Y_piers, Z_piers, (num_spandrels+1:num_spandrels+num_piers)', is_pier];

        % Add restraints
        patterns = {'ry+rz', 'rx+ry'};
        % Convert data to a cell array
        piers = num2cell(piers);
        % Initialize the counter for the pattern index
        pattern_index = 1;
        % Loop through the rows of data
        for i = 1:size(piers, 1)
            % Get the pattern for the current row
            pattern = patterns{pattern_index};
            % Assign the pattern to the 6th column of the current row
            piers{i, 6} = pattern;
            if mod(i, 4) == 0
                % Update the pattern index for every 4 rows
                pattern_index = mod(pattern_index, numel(patterns)) + 1;
            end
        end
        % Sort the piers combined matrix based on Z first, then Y, and finally on X within each Z-Y group
        sorted_piers = sortrows(piers, [3, 2, 1]);

        % Create a column vector to indicate whether each element is a piers (1) or not (0)
        is_node = ones(num_nodes, 1).*2; % Assuming all elements are not piers initially
        % Create a combined matrix of XPiers, YPiers, ZPiers, their indices, and the piers indicator
        nodes= [X_nodes, Y_nodes, Z_nodes, (num_spandrels+num_piers+1:num_elements)', is_node];

        % Add restraints
        patterns = {'ry+rz', 'rx+ry'};
        % Convert data to a cell array
        nodes = num2cell(nodes);
        % Initialize the counter for the pattern index
        pattern_index = 1;
        % Loop through the rows of data
        for i = 1:size(nodes, 1)
            % Get the pattern for the current row
            pattern = patterns{pattern_index};
            % Assign the pattern to the 6th column of the current row
            nodes{i, 6} = pattern;
            if mod(i, 4) == 0
                % Update the pattern index for every 4 rows
                pattern_index = mod(pattern_index, numel(patterns)) + 1;
            end
        end

        % Sort the piers combined matrix based on Z first, then Y, and finally on X within each Z-Y group
        sorted_nodes = sortrows(nodes, [3, 2, 1]);

        % Concatenate the combined matrices for spandrels and piers
        combined = [spandrels; piers; sorted_nodes];

        % Sort the combined matrix based on Z first, then Y, and finally on X within each Z-Y group
        sorted_combined = [sorted_spandrels; sorted_piers; sorted_nodes];

        % Extract sorted X, Y, Z coordinates
        X_sorted = cell2mat(sorted_combined(:, 1));
        Y_sorted = cell2mat(sorted_combined(:, 2));
        Z_sorted = cell2mat(sorted_combined(:, 3));

        count=2;
        Xidx=X_sorted(1);
        for i=2:length(X_sorted)
            if X_sorted(i)~=X_sorted(i-1)
                Xidx(count)=X_sorted(i);
                count=count+1;
            end
        end

        count=2;
        Zidx=Z_sorted(1);
        for i=2:length(Z_sorted)
            if Z_sorted(i)~=Z_sorted(i-1)
                Zidx(count)=Z_sorted(i);
                count=count+1;
            end
        end

        % Generate node names based on sorted indices
        node_name = cell(size(X_sorted));
        x_count=1;
        z_count=1;
        type=0; %spandrel(0), pier(1) or node(2)
        z_level_prev=1;
        x_level_prev=1;

        for i = 1:numel(X_sorted)

            % Extracting the indices of z for naming convention
            if Zidx(z_count) == Z_sorted(i)
                z_level = z_level_prev; % Get the z level

            elseif Zidx(z_count) > Z_sorted(i) %indicates change of element type
                z_count=z_count+1;
                z_level_prev=1;
                z_level = z_level_prev; % Get the z level
                type=type+1;
            else % change of level within same element type
                z_count=z_count+1;
                z_level = z_level_prev+1;
                z_level_prev=z_level;
                x_level_prev=0; %indicates change in z value
            end

            % Extracting the indices of x for naming convention
            if Xidx(x_count) == X_sorted(i)
                x_level = x_level_prev; % Get the x level
            else
                x_count=x_count+1;
                x_level = x_level_prev+1;
                x_level_prev=x_level;
            end

            % Generate the node name based on the pattern
            node_name{i} = sprintf('%d%d%d%d%d%d%d', z_level, 0, type, 0, i);
        end

        % Convert node_name to a column vector
        node_name = node_name(:);

        % Convert numeric arrays to cell arrays
        X_sorted_cell = num2cell(X_sorted);
        Y_sorted_cell = num2cell(Y_sorted);
        Z_sorted_cell = num2cell(Z_sorted);
        index_cell= num2cell(sorted_combined(:,4));
        s_or_p_cell = num2cell(sorted_combined(:, 5));

        % Combine sorted data with node names
        sorted_data_with_names = [X_sorted_cell, Y_sorted_cell, Z_sorted_cell, index_cell, s_or_p_cell, node_name];

        % Match node names with original combine
        % d matrix indices
        combined_with_names = cell(size(combined, 1), size(combined, 2) + 1); % Initialize combined_with_names
        combined_with_names(:, 1:end-1) = num2cell(combined); % Copy X, Y, Z, indices to combined_with_names

        % Match node names with corresponding indices
        for i = 1:size(combined, 1)
            idx = find(cell2mat(sorted_combined(:, 4)) == i, 1); % Get the index in the sorted_combined matrix
            combined_with_names{i, end} = node_name{idx}; % Assign the corresponding node name
        end



        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ELEMENT CONNECTIVITY %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Define the number of nodes per element
        nodes_per_element = 8;
        % Determine the number of elements
        num_elements = size(combined_with_names, 1) / nodes_per_element;
        % Initialize the matrix to store elements and their node names
        elements_with_nodes = cell(num_elements, 3);

        % Extract node names for each element
        for i = 1:num_elements
            % Determine the indices for the current element
            start_idx = (i - 1) * nodes_per_element + 1;
            end_idx = start_idx + nodes_per_element - 1;

            % Extract node names for the current element
            node_names = combined_with_names(start_idx:end_idx, end);

            % Assign element name
            elem_name = sprintf('00%d', i); % should be digits only

            % Store node names for the current element
            elements_with_nodes{i} = node_names';

            % Store element name and node names for the current element
            elements_with_nodes{i, 1} = elem_name;
            elements_with_nodes{i, 2} = 'macro1';
            elements_with_nodes{i, 3} = node_names';
        end

        % Display the elements with their node names
        for i = 1:num_elements
            fprintf('Element %d (%s): %s\n', i, elements_with_nodes{i, 1}, strjoin(elements_with_nodes{i, 3}, ' '));
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% INITIAL LOADS %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % check if all values of array are the same, in which case this
            % is a W_Y case, therefore the idx =2
        if numel(unique(cellfun(@(x) x, sorted_data_with_names(:, 1)))) == 1
            idx=2;
        else
            idx=1;
        end

        init_load_sorted_data= sortrows(sorted_data_with_names, [idx, 3]);

        % Find unique values in the idx column
        unique_first_column = unique(cell2mat(init_load_sorted_data(:, idx)));

        % Initialize the init_load_nod vector
        init_load_nod = cell(length(unique_first_column), 1);

        % Iterate over unique values in the first column
        for i = 1:length(unique_first_column)
            % Find rows corresponding to the current unique value in the idx column
            rows = find(cell2mat(init_load_sorted_data(:, idx)) == unique_first_column(i));

            % Extract the subset of data for the current group (of same x node values)
            group_data = init_load_sorted_data(rows, :);

            % Find the maximum value in the third column for the current group
            [~, index] = max(cell2mat(group_data(:, 3)));

            % Store the corresponding value from the sixth column in init_load_nod
            init_load_nod{i} = group_data{index, 6};
        end

        % Initialize the vector 'direction'
        direction = cell(length(init_load_nod), 1);

        % Fill the vector of strings 'y'  
        for i = 1:length(init_load_nod)
            direction{i} = 'y';
        end

        % Initialize the vector 'type'
        type = cell(length(init_load_nod), 1);

        % Fill the vector of strings 'f'
        for i = 1:length(init_load_nod)
            type{i} = 'f';
        end

        % Initialize the vector of doubles with zeros
        vector_length = length(init_load_nod);
        value = zeros(vector_length, 1);

        % Fill the vector with the value -25000
        value(:) = -25000;

        %% Append +z restraint to top nodes

        for i = 1:length(combined_with_names(:, 7))
            % Check if the current element of 7th col is in init_load_nod
            if any(strcmp(combined_with_names{i, 7}, init_load_nod))==1
                % If it is, append '+z' to the corresponding element in the 6th column
                combined_with_names{i, 6} = append(combined_with_names{i, 6}, '+z');
            end
        end



        %%%%%%%%%%%%%%%%%%%%
        %%%% SAVE FILES %%%%
        %%%%%%%%%%%%%%%%%%%%
        save_file(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls/sorted_coordinates.txt', X_sorted, Y_sorted, Z_sorted, node_name, 'sorted_coordinates');
        save_file(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls/elem_coordinates.txt', num_elements, elements_with_nodes, 'elem_coordinates');
        save_file(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls/restraints.txt', combined_with_names, 'restraints');
        save_file(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls/initial_loads.txt', init_load_nod, direction, type, value, 'initial_loads');

%                 % Save X_sorted, Y_sorted, Z_sorted, and node_name to a text file on the desktop
%                 output_file = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/sorted_coordinates.txt';
%                 fid = fopen(output_file, 'w');
%                 if fid == -1
%                     error('Could not open file for writing: %s', output_file);
%                 end
%         
%                 fprintf(fid, 'nod.name\tx\ty\tz\n');
%                 for i = 1:numel(X_sorted)
%                     fprintf(fid, '%s\t%f\t%f\t%f\n', node_name{i}, X_sorted(i), Y_sorted(i), Z_sorted(i));
%                 end
%         
%         
%                 % Save elements_with_nodes to a text file on the desktop
%                 output_file = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/elem_coordinates.txt';
%                 fid = fopen(output_file, 'w');
%                 if fid == -1
%                     error('Could not open file for writing: %s', output_file);
%                 end
%         
%                 fprintf(fid, 'elm.name\tgrp.name \tnod.name\n');
%                 for i = 1:num_elements
%                     % Get the element name and node names of the current element
%                     elem_name = elements_with_nodes{i, 1};
%                     macro_name = elements_with_nodes{i, 2};
%                     node_names = elements_with_nodes{i, 3};
%         
%                     % Concatenate node names of the current element with spaces
%                     node_names_str = strjoin(node_names, ' ');
%         
%                     % Print the element name, 'macro1', and concatenated node names
%                     fprintf(fid, '%s\t%s\t%s\n', elem_name, macro_name, node_names_str);
%                 end
%         
%                 % Save node names and restraints to a text file on the desktop
%                 output_file = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/restraints.txt';
%                 fid = fopen(output_file, 'w');
%                 if fid == -1
%                     error('Could not open file for writing: %s', output_file);
%                 end
%         
%                 fprintf(fid, 'nod.name\tdirection\n');
%                 for i = 1:numel(combined_with_names(:,1))
%                     node_names=combined_with_names{i, 7};
%                     restraints=combined_with_names{i, 6};
%                     fprintf(fid, '%s\t%s\n', string(node_names), string(restraints));
%                 end
%         
%         
%                 % Save initial loads to a text file
%                 output_file = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/initial_loads.txt';
%                 fid = fopen(output_file, 'w');
%                 if fid == -1
%                     error('Could not open file for writing: %s', output_file);
%                 end
%         
%                 fprintf(fid, 'nod.name\tdirection\ttype\tvalue\n');
%                 for i = 1:numel(init_load_nod)
%                     fprintf(fid, '%s\t%s\t%s\t%d\n', init_load_nod{i}, direction{i}, type{i}, value(i));
%                 end
%         
%         
%                 fclose(fid);
                 successOutput = true;

    end
end

end


