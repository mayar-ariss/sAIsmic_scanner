% Mayar Ariss 10 May 2024

function [successOutput, node_name] = wallAdaptic_disc(XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face)

addpath('/Users/mayar/Desktop/ICL/4Y/Thesis/Codes/FEM_buildings/src/utils/hex_and_rgb_v1.1.1'); %careful, / for macOS and \ for Windows

%FILE NAME
fn='_BUILDING2';

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
        X_nodes = round(XNodes_filtered(:, zero_indices), 1);
        Y_nodes = round(YNodes_filtered(:, zero_indices), 1);
        Z_nodes = round(ZNodes_filtered(:, zero_indices), 1);

        % Convert to vectors
        X_nodes = X_nodes(:);
        Y_nodes = Y_nodes(:);
        Z_nodes = Z_nodes(:);


        %% SPANDRELS

        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XSpandrels_filtered, 4) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YSpandrels_filtered, 4) == values);
        end

        % Select X, Y, and Z for only front frame
        X_spandrels = round(XSpandrels_filtered(:, zero_indices), 1);
        Y_spandrels = round(YSpandrels_filtered(:, zero_indices), 1);
        Z_spandrels = round(ZSpandrels_filtered(:, zero_indices), 1);

        % Reshape X, Y, and Z matrices to column vectors
        X_spandrels = X_spandrels(:);
        Y_spandrels = Y_spandrels(:);
        Z_spandrels = Z_spandrels(:);


        %% PIERS

        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XPiers_filtered, 4) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YPiers_filtered, 4) == values);
        end


        X_piers = round(XPiers_filtered(:, zero_indices), 1);
        Y_piers = round(YPiers_filtered(:, zero_indices), 1);
        Z_piers = round(ZPiers_filtered(:, zero_indices), 1);

        % Reshape XPiers, YPiers, and ZPiers matrices to column vectors
        X_piers = X_piers(:);
        Y_piers = Y_piers(:);
        Z_piers = Z_piers(:);


        %% BLOCKS

        % All elements
        X_elems= [X_nodes; X_spandrels; X_piers];
        Y_elems= [Y_nodes; Y_spandrels; Y_piers];
        Z_elems= [Z_nodes; Z_spandrels; Z_piers];


        % DATA CLEANING
        X_blocks= X_elems';
        Y_blocks= Y_elems';
        Z_blocks= Z_elems';


        % Rearrange every four-element cluster
        for i = 1:8:numel(X_blocks)
            % Extract the current four-element cluster
            cluster_x = X_blocks(i:i+3);
            cluster_y = Y_blocks(i:i+3);
            cluster_z = Z_blocks(i:i+3);

            % Sort the cluster
            sorted_cluster_x = sort(cluster_x);
            sorted_cluster_y = sort(cluster_y);
            sorted_cluster_z = sort(cluster_z);

            % Rearrange the cluster as a,b,b,a where a<b
            X_blocks(i:i+3) = sorted_cluster_x([2, 3, 3, 2]);
            X_blocks(i+4:i+7) = sorted_cluster_x([3, 2, 2, 3]);

            Y_blocks(i:i+3) = sorted_cluster_y([2, 3, 3, 2]);
            Y_blocks(i+4:i+7) = sorted_cluster_y([3, 2, 2, 3]);

            Z_blocks(i:i+3) = sorted_cluster_z([2, 2, 3, 3]);
            Z_blocks(i+4:i+7) = sorted_cluster_z([2, 2, 3, 3]);

        end

        X_elems=X_blocks';
        Y_elems=Y_blocks';
        Z_elems=Z_blocks';

        % Number of block elements (/4 to account for number of sides, includes cohesive interfaces)
        num_blocks = size(X_elems, 1)/4;
        j=0;

        % Initialize block variable to store coordinates
        blocks = zeros(4, 3, 1, num_blocks); % 4 rows for x and y coordinates, 2 columns for x and y

        % Iterate over each block element
        for i = 1:4:length(X_elems)

            %keep count of n° of blocks
            j=j+1;

            % Extract x and y coordinates for the current block element
            x_coords = X_elems(i:i+3);
            y_coords = Y_elems(i:i+3);
            z_coords = Z_elems(i:i+3);

            % Store coordinates in the blocks variable
            if numel(unique(X_elems(:))) == 1
                % 2 walls where x_coords is constant
                blocks(:,:,:,j) = [y_coords, z_coords, x_coords];
            elseif numel(unique(Y_elems(:))) == 1
                % 2 walls where y_coords is constant
                blocks(:,:,:,j) = [x_coords, z_coords, y_coords];
            end


        end


        %% HORIZONTAL DISCRETISATION

        % Initialize new_blocks initial
        new_blocks_init = zeros([4 3 1 0]);

        %initialise new_blocks_final
        new_blocks_final=blocks;

        while size(new_blocks_init, 4) ~= size(new_blocks_final, 4)

            %initialise added_blocks
            added_blocks=0;

            %number of initial blocks
            new_blocks_init=new_blocks_final;

            %number of blocks
            num_blocks= size(new_blocks_init, 4);

            % Calling this before the loop to initialize.
            finishNow = false;  

            for i=1:num_blocks

                % Select the coordinates of the first block
                my_block = new_blocks_init(:,:,:,i);

                % Given rectangle coordinates
                x = my_block(:, 1);
                y = my_block(:, 2);

                % Initialize variables to store side coordinates
                side1 = zeros(2, 2);
                side2 = zeros(2, 2);
                side3 = zeros(2, 2);
                side4 = zeros(2, 2);

                % Loop through each side of the rectangle
                for side = 1:4
                    % Determine the indices of the vertices forming the current side
                    start_idx = side;
                    end_idx = mod(side, 4) + 1;

                    % Extract coordinates of the vertices forming the current side
                    side_x = [x(start_idx), x(end_idx)];
                    side_y = [y(start_idx), y(end_idx)];

                    % Store side coordinates in corresponding variables
                    switch side
                        case 1
                            side1 = [side_x', side_y'];
                        case 2
                            side2 = [side_x', side_y'];
                        case 3
                            side3 = [side_x', side_y'];
                        case 4
                            side4 = [side_x', side_y'];
                    end
                end

                % Check if any of the vertices of the current block are on any side of other blocks
                for j = 1:num_blocks

                    other_block = new_blocks_init(:,:,:,j);

                    % Check if any vertex of any other block lies on any side of the current block
                    for k = 1:4
                        vertex = other_block(k,:);
                        for side = 1:4
                            my_side = my_block([side, mod(side, 4) + 1], :);

                            if vertex(1,1)==my_side(1,1) && vertex(1,1)== my_side(2,1) ...
                                    && min(my_side(:,2)) < vertex(1,2) && vertex(1,2) < max(my_side(:,2))
                                disp(['FACE ', num2str(values),' Vertex ', num2str(k), ' of my block ', num2str(i), ' lies on side ', num2str(side), ' of other block ', num2str(j)]);

                                % new block
                                new_block1= [my_block(:, 1), [min(my_block(:,2)), min(my_block(:,2)), vertex(1,2), vertex(1,2)]', my_block(:, 3)];
                                % new cohesive interface
                                new_block2= [my_block(:, 1), [min(my_block(:,2)), min(my_block(:,2)), vertex(1,2), vertex(1,2)]', my_block(:, 3)];
                                % new block
                                new_block3= [my_block(:, 1), [vertex(1,2), vertex(1,2), max(my_block(:,2)), max(my_block(:,2))]', my_block(:, 3)];
                                % new cohesive interface
                                new_block4= [my_block(:, 1), [vertex(1,2), vertex(1,2), max(my_block(:,2)), max(my_block(:,2))]', my_block(:, 3)];

                                % delete cohesive block
                                new_blocks_final(:,:,:,i+added_blocks)=[];

                                % Replace the deleted block with new_block1
                                new_blocks_final(:,:,:,i+added_blocks) = new_block1;

                                % Determine the index for inserting new_block2
                                new_index = i + added_blocks;
                                % Insert new_block2 at the determined index
                                new_blocks_final = cat(4, new_blocks_final(:,:,:,1:new_index-1), new_block2, new_blocks_final(:,:,:,new_index:end));

                                % Determine the index for inserting new_block3
                                new_index = new_index + 1 + 1;
                                % Insert new_block3 at the determined index
                                new_blocks_final = cat(4, new_blocks_final(:,:,:,1:new_index-1), new_block3, new_blocks_final(:,:,:,new_index:end));

                                % Determine the index for inserting new_block4
                                new_index = new_index + 1;
                                % Insert new_block4 at the determined index
                                new_blocks_final = cat(4, new_blocks_final(:,:,:,1:new_index-1), new_block4, new_blocks_final(:,:,:,new_index:end));

                                % keep track of number of blocks added (2, 3 and 4)
                                added_blocks=added_blocks+2;

                                finishNow=true;

                            elseif vertex(1,2)==my_side(1,2) && vertex(1,2)==my_side(2,2) ...
                                    && min(my_side(:,1)) < vertex(1,1) && vertex(1,1) < max(my_side(:,1))
                                disp(['vert ', 'FACE ', num2str(values),' Vertex ', num2str(k), ' of my block ', num2str(i), ' lies on side ', num2str(side), ' of other block ', num2str(j)]);

                                % new block
                                new_block1= [[min(my_block(:,1)), vertex(1,1), vertex(1,1), min(my_block(:,1))]', my_block(:, 2),  my_block(:, 3)];
                                % new cohesive interface
                                new_block2= [[vertex(1,1), min(my_block(:,1)), min(my_block(:,1)), vertex(1,1)]', my_block(:, 2), my_block(:, 3)];
                                % new block
                                new_block3= [[vertex(1,1), max(my_block(:,1)), max(my_block(:,1)), vertex(1,1)]', my_block(:, 2), my_block(:, 3)];
                                % new cohesive interface
                                new_block4= [[max(my_block(:,1)), vertex(1,1), vertex(1,1), max(my_block(:,1))]', my_block(:, 2), my_block(:, 3)];

                                % delete cohesive block
                                new_blocks_final(:,:,:,i+added_blocks)=[];

                                % Replace the deleted block with new_block1
                                new_blocks_final(:,:,:,i+added_blocks) = new_block1;

                                % Determine the index for inserting new_block2
                                new_index = i + added_blocks;
                                % Insert new_block2 at the determined index
                                new_blocks_final = cat(4, new_blocks_final(:,:,:,1:new_index-1), new_block2, new_blocks_final(:,:,:,new_index:end));

                                % Determine the index for inserting new_block3
                                new_index = new_index + 1 + 1;
                                % Insert new_block3 at the determined index
                                new_blocks_final = cat(4, new_blocks_final(:,:,:,1:new_index-1), new_block3, new_blocks_final(:,:,:,new_index:end));

                                % Determine the index for inserting new_block4
                                new_index = new_index + 1;
                                % Insert new_block4 at the determined index
                                new_blocks_final = cat(4, new_blocks_final(:,:,:,1:new_index-1), new_block4, new_blocks_final(:,:,:,new_index:end));

                                % keep track of number of blocks added (2, 3 and 4)
                                added_blocks=added_blocks+2;

                                finishNow=true;

                            end

                            if finishNow
                                break
                            end

                        end

                        if finishNow
                            break
                        end

                    end

                    if finishNow
                        break
                    end

                end

                if finishNow
                    break
                end

            end

            debugging=1;

        end



        %% ADD RESTRAINTS

        X_blocks= reshape(squeeze(new_blocks_final(:, 1, :)), [], 1);
        Y_blocks= reshape(squeeze(new_blocks_final(:, 2, :)), [], 1);
        Z_blocks= reshape(squeeze(new_blocks_final(:, 3, :)), [], 1);


        % Rearrange every four-element cluster
        for i = 1:8:numel(X_blocks)
            % Extract the current four-element cluster
            cluster_x = X_blocks(i:i+3);
            cluster_y = Y_blocks(i:i+3);
            cluster_z = Z_blocks(i:i+3);

            % Sort the cluster
            sorted_cluster_x = sort(cluster_x);
            sorted_cluster_y = sort(cluster_y);
            sorted_cluster_z = sort(cluster_z);

            % Rearrange the cluster as a,b,b,a where a<b
            X_blocks(i:i+3) = sorted_cluster_x([2, 3, 3, 3]);
            X_blocks(i+4:i+7) = sorted_cluster_x([3, 2, 2, 2]);

            Y_blocks(i:i+3) = sorted_cluster_y([2, 2, 2, 3]);
            Y_blocks(i+4:i+7) = sorted_cluster_y([3, 3, 3, 2]);

            Z_blocks(i:i+3) = sorted_cluster_z([2, 3, 3, 3]);
            Z_blocks(i+4:i+7) = sorted_cluster_z([3, 2, 2, 2]);

        end

        % Determine the number of total elements
        num_elements=numel(X_blocks);

        % Create a combined matrix of X, Y, Z, and their indices
        block_mat = [X_blocks, Y_blocks, Z_blocks, (1:num_elements)'];

        patterns = {'ry+rz', 'rx+ry'};
        % Convert data to a cell array
        block_mat = num2cell(block_mat);
        % Initialize the counter for the pattern index
        pattern_index = 1;
        % Loop through the rows of data
        for i = 1:size(block_mat, 1)
            % Get the pattern for the current row
            pattern = patterns{pattern_index};
            % Assign the pattern to the 6th column of the current row
            block_mat{i, 5} = pattern;
            if mod(i, 4) == 0
                % Update the pattern index for every 4 rows (alternative between cohesive interfaces and blocks)
                pattern_index = mod(pattern_index, numel(patterns)) + 1;
            end
        end
        % Sort the blocks matrix based on Z first, then Y, and finally on X within each Z-Y group
        sorted_blocks = sortrows(block_mat, [3, 2, 1]);

        % Extract sorted X, Y, Z coordinates
        X_sorted = cell2mat(sorted_blocks(:, 1));
        Y_sorted = cell2mat(sorted_blocks(:, 2));
        Z_sorted = cell2mat(sorted_blocks(:, 3));

        count=2;
        Xidx=X_sorted(1);
        for i=2:length(X_sorted)
            if X_sorted(i)~=X_sorted(i-1)
                Xidx(count)=X_sorted(i);
                count=count+1;
            end
        end

        count=2;
        Yidx=Y_sorted(1);
        for i=2:length(Y_sorted)
            if Y_sorted(i)~=Y_sorted(i-1)
                Yidx(count)=Y_sorted(i);
                count=count+1;
            end
        end

        % Generate node names based on sorted indices
        node_name = cell(size(X_sorted));
        x_count=1;
        y_count=1;
        y_level_prev=1;
        x_level_prev=1;

        for i = 1:numel(X_sorted)

            % Extracting the indices of z for naming convention
            if Yidx(y_count) == Y_sorted(i)
                y_level = y_level_prev; % Get the z level

            else % change of level within same element type
                y_count=y_count+1;
                y_level = y_level_prev+1;
                y_level_prev=y_level;
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
            node_name{i} = sprintf('%d%d%d%d%d%d%d', y_level, 0, 0, i);
        end

        % Convert node_name to a column vector
        node_name = node_name(:);

        % Convert numeric arrays to cell arrays
        X_sorted_cell = num2cell(X_sorted);
        Y_sorted_cell = num2cell(Y_sorted);
        Z_sorted_cell = num2cell(Z_sorted);
        index_cell= num2cell(sorted_blocks(:,4));

        % Combine sorted data with node names
        sorted_data_with_names = [X_sorted_cell, Y_sorted_cell, Z_sorted_cell, index_cell, node_name];

        % Match node names with original combine
        % d matrix indices
        blocks_with_names = cell(size(sorted_blocks, 1), size(sorted_blocks, 2) + 1); % Initialize combined_with_names
        blocks_with_names(:, 1:end-1) = num2cell(sorted_blocks); % Copy X, Y, Z, indices to combined_with_names

        % Match node names with corresponding indices
        for i = 1:size(sorted_blocks, 1)
            idx = find(cell2mat(sorted_blocks(:, 4)) == i, 1); % Get the index in the sorted_combined matrix
            blocks_with_names{i, end} = node_name{idx}; % Assign the corresponding node name
        end


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ELEMENT CONNECTIVITY %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Define the number of nodes per element
        nodes_per_element = 8;
        % Determine the number of elements
        num_elements = size(blocks_with_names, 1) / nodes_per_element;
        % Initialize the matrix to store elements and their node names
        elements_with_nodes = cell(num_elements, 3);

        % Extract node names for each element
        for i = 1:num_elements
            % Determine the indices for the current element
            start_idx = (i - 1) * nodes_per_element + 1;
            end_idx = start_idx + nodes_per_element - 1;

            % Extract node names for the current element
            node_names = blocks_with_names(start_idx:end_idx, end);

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

            % Find the maximum value in the second column y for the current group
            [~, index] = max(cell2mat(group_data(:, 2)));

            % Store the corresponding value from the sixth column in init_load_nod
            init_load_nod{i} = group_data{index, 5};
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

        for i = 1:length(blocks_with_names(:, end))
            % Check if the current element of 7th col is in init_load_nod
            if any(strcmp(blocks_with_names{i, 6}, init_load_nod))
                % If it is, append '+z' to the corresponding element in the 6th column
                blocks_with_names{i, 5} = append(blocks_with_names{i, 5}, '+z');
            end
        end



        %%%%%%%%%%%%%%%%%%%%
        %%%% SAVE FILES %%%%
        %%%%%%%%%%%%%%%%%%%%
        save_file_disc(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls_disc/sorted_coordinates.txt', X_sorted, Y_sorted, Z_sorted, node_name, 'sorted_coordinates');
        save_file_disc(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls_disc/elem_coordinates.txt', num_elements, elements_with_nodes, 'elem_coordinates');
        save_file_disc(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls_disc/restraints.txt', blocks_with_names, 'restraints');
        save_file_disc(prefix, '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/all_walls_disc/initial_loads.txt', init_load_nod, direction, type, value, 'initial_loads');





        %reset added_blocks
        added_blocks=0;

    end

end
successOutput = true;

end
