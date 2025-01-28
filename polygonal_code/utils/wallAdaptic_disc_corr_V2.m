% Mayar Ariss 10 May 2024

function [successOutput] = wallAdaptic_disc_corr_V2(XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face)

addpath('/Users/mayar/Desktop/ICL/4Y/Thesis/Codes/FEM_buildings/src/utils/hex_and_rgb_v1.1.1'); %careful, / for macOS and \ for Windows

%FILE NAME
fn='_BUILDING2_corr';
save_file_path = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/building_2_disc/sorted_coordinates.txt';

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

        patterns = {'ry+rz', 'rx+rz'};
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
        blocks_with_names = cell(size(sorted_data_with_names, 1), size(sorted_data_with_names, 2) + 1); % Initialize combined_with_names
        blocks_with_names(:, 1:end-1) = block_mat; % Copy X, Y, Z, indices to combined_with_names

        % Match node names with corresponding indices
        for i = 1:size(sorted_blocks, 1)
            idx = find(cell2mat(sorted_blocks(:, 4)) == i, 1); % Get the index in the sorted_blocks matrix
            blocks_with_names{i, end} = node_name{idx}; % Assign the corresponding node name
        end


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ELEMENT CONNECTIVITY %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Final Mesh

        meshmat=blocks_with_names;

        % Extract the unique coordinates
        [unique_coords, ~, idx] = unique(cell2mat(meshmat(:, 1:2)), 'rows', 'stable');

        % Find groups of nodes sharing the same coordinates
        groups = accumarray(idx, 1:numel(idx), [], @(x) {sort(x)});

        % Iterate through groups
        for i = 1:numel(groups)

            group_indices = groups{i};

            if numel(group_indices) == 4 % Check if group has four nodes

                % Variables to store nodes based on fifth column condition
                rx_rz_gpis = {};
                ry_rz_gpis = {};

                for j = 1:numel(group_indices)
                    if strcmp(meshmat{group_indices(j), 5}, 'rx+rz')
                        rx_rz_gpis = [rx_rz_gpis; group_indices(j)];
                    elseif strcmp(meshmat{group_indices(j), 5}, 'ry+rz')
                        ry_rz_gpis = [ry_rz_gpis; group_indices(j)];
                    end
                end

                if numel(rx_rz_gpis) == 3
                    meshmat{rx_rz_gpis{2}, 6} = [meshmat{rx_rz_gpis{1}, 6}];
                end
                if numel(ry_rz_gpis) == 3
                    meshmat{ry_rz_gpis{2}, 6} = [meshmat{ry_rz_gpis{1}, 6}];
                end

            elseif numel(groups{i}) == 6 % Check if the group has six nodes

                rx_rz_gpis = {};
                ry_rz_gpis = {};

                for j = 1:numel(group_indices)
                    if strcmp(meshmat{group_indices(j), 5}, 'rx+rz')
                        rx_rz_gpis = [rx_rz_gpis; group_indices(j)];
                    elseif strcmp(meshmat{group_indices(j), 5}, 'ry+rz')
                        ry_rz_gpis = [ry_rz_gpis; group_indices(j)];
                    end
                end

                if numel(rx_rz_gpis) == 3
                    meshmat{rx_rz_gpis{2}, 6} = [meshmat{rx_rz_gpis{1}, 6}];
                end

                if numel(ry_rz_gpis) == 3
                    meshmat{ry_rz_gpis{2}, 6} = [meshmat{ry_rz_gpis{1}, 6}];
                    %meshmat{rx_rz_gpis{3}, 6} = [meshmat{ry_rz_gpis{3}, 6}];
                end
                %%%
                if numel(rx_rz_gpis) == 4
                    meshmat{rx_rz_gpis{3}, 6} = [meshmat{rx_rz_gpis{1}, 6}];
                    %meshmat{rx_rz_gpis{4}, 6} = [meshmat{rx_rz_gpis{2}, 6}];
                end
                if numel(ry_rz_gpis) == 2
                    meshmat{ry_rz_gpis{2}, 6} = [meshmat{ry_rz_gpis{1}, 6}];
                end
                %%%
                if numel(ry_rz_gpis) == 4
                    meshmat{ry_rz_gpis{3}, 6} = [meshmat{ry_rz_gpis{1}, 6}];
                    %meshmat{ry_rz_gpis{4}, 6} = [meshmat{ry_rz_gpis{2}, 6}];
                end
                if numel(rx_rz_gpis) == 2
                    meshmat{rx_rz_gpis{2}, 6} = [meshmat{rx_rz_gpis{1}, 6}];
                end
                %%%
          
            elseif numel(groups{i}) == 8 % Check if the group has eight nodes


                rx_rz_gpis = {};
                ry_rz_gpis = {};

                for j = 1:numel(group_indices)
                    if strcmp(meshmat{group_indices(j), 5}, 'rx+rz')
                        rx_rz_gpis = [rx_rz_gpis; group_indices(j)];
                    elseif strcmp(meshmat{group_indices(j), 5}, 'ry+rz')
                        ry_rz_gpis = [ry_rz_gpis; group_indices(j)];
                    end
                end

                if numel(rx_rz_gpis) == 4
                    meshmat{rx_rz_gpis{3}, 6} = [meshmat{rx_rz_gpis{1}, 6}];
                    meshmat{rx_rz_gpis{4}, 6} = [meshmat{rx_rz_gpis{2}, 6}];
                end
                if numel(ry_rz_gpis) == 4
                    meshmat{ry_rz_gpis{3}, 6} = [meshmat{ry_rz_gpis{1}, 6}];
                    meshmat{ry_rz_gpis{4}, 6} = [meshmat{ry_rz_gpis{2}, 6}];
                end

            end
        end



        % Define the number of nodes per element
        nodes_per_element = 8;
        % Determine the number of elements
        num_elements = size(meshmat, 1) / nodes_per_element;
        % Initialize the matrix to store elements and their node names
        elements_with_nodes = cell(num_elements, 3);



        %         % detect top elements which will be defined as macro2
        %         % check if all values of array are the same, in which case this
        %         % is a W_Y case, therefore the idx =2
        %         if numel(unique(cellfun(@(x) x, meshmat(:, 1)))) == 1
        %             idx=2;
        %         else
        %             idx=1;
        %         end
        %
        %         init_restraint_sorted_data= sortrows(meshmat, [idx, 3]);
        %
        %         % Find unique values in the idx column
        %         unique_first_column = unique(cell2mat(init_restraint_sorted_data(:, idx)));
        %
        %         % Initialize the init_load_nod vector
        %         init_rigid_nod = {};
        %
        %         % Iterate over unique values in the first column
        %         for i = 1:length(unique_first_column)
        %             % Find rows corresponding to the current unique value in the idx column
        %             rows = find(cell2mat(init_restraint_sorted_data(:, idx)) == unique_first_column(i));
        %
        %             % Extract the subset of data for the current group (of same x node values)
        %             group_data = init_restraint_sorted_data(rows, :);
        %
        %             % Find the maximum value(s) in the second column y for the current group
        %             max_values = max(cell2mat(group_data(:, 2)));
        %
        %             % Find the rows where the second column matches the maximum value(s)
        %             max_rows = ismember(cell2mat(group_data(:, 2)), max_values);
        %
        %             % Store the corresponding value from the sixth column in init_load_nod
        %             init_rigid_nod = [init_rigid_nod; group_data(max_rows, 6)];
        %         end





        % Extract node names for each element
        for i = 1:num_elements
            % Determine the indices for the current element
            start_idx = (i - 1) * nodes_per_element + 1;
            end_idx = start_idx + nodes_per_element - 1;

            % Extract node names for the current element
            node_names = meshmat(start_idx:end_idx, end);

            % Assign element name
            elem_name = sprintf('100%d', i); % should be digits only

            % Store node names for the current element
            elements_with_nodes{i} = node_names';

            % Store element name and node names for the current element
            elements_with_nodes{i, 1} = elem_name;

            % Check if any string in node_names starts with '1'
            if any(cellfun(@(x) startsWith(x, '1'), node_names)) %|| ...
                % any(ismember(node_names, init_rigid_nod))
                % If any string in node_names starts with '1', assign 'macro2'
                elements_with_nodes{i, 2} = 'macro2';
            else
                % If no string in node_names starts with '1', assign 'macro1'
                elements_with_nodes{i, 2} = 'macro1';
            end

            elements_with_nodes{i, 3} = node_names';
        end

        % Display the elements with their node names
        for i = 1:num_elements
            fprintf('Element %d (%s): %s\n', i, elements_with_nodes{i, 1}, strjoin(elements_with_nodes{i, 3}, ' '));
        end


        %% redefine mesh nodal coordinates

        %delete duplicate rows to account for shared sides
        [~,uidx] = unique(meshmat(:,6),'stable');
        meshmat_final = meshmat(uidx,:);

        %reset new index
        idx_final=1:size(meshmat_final, 1);
        meshmat_final(:, 4)=num2cell(idx_final);

        %sort the matrix again
        meshmat_final=sortrows(meshmat_final, [2 1]);

        %the vectors used for ADAPTIC file
        X_final=cell2mat(meshmat_final(:,1));
        Y_final=cell2mat(meshmat_final(:,2));
        Z_final=cell2mat(meshmat_final(:,3));
        node_names_final=meshmat_final(:,6);


        % Change the restraints according to final mesh
        % Extract the unique coordinates
        [unique_coords, ~, idx] = unique(cell2mat(meshmat_final(:, 1:2)), 'rows', 'stable');

        % groups of nodes sharing the same coordinates
        groups_final = accumarray(idx, 1:numel(idx), [], @(x) {sort(x)});

        % Iterate through groups
        for i = 1:numel(groups_final)
            group_indices = groups_final{i};
            if numel(group_indices) == 2 % Check if group has two nodes
                % Append node restraints as specified
                meshmat_final{group_indices(1), 5} = 'ry+rz';
                meshmat_final{group_indices(2), 5} = 'rx+rz';

            elseif numel(groups_final{i}) == 3 % Check if the group has six nodes
                % Replace node names for groups of six
                meshmat_final{groups_final{i}(1), 5} = 'ry+rz'; % Replace first with third
                meshmat_final{groups_final{i}(2), 5} = 'rx+rz'; % Replace second with fourth
                meshmat_final{groups_final{i}(3), 5} = 'ry+rz'; % Replace first with third

            elseif numel(groups_final{i}) == 4 % Check if the group has eight nodes
                % Replace node names for groups of eight
                meshmat_final{groups_final{i}(1), 5} = 'ry+rz'; % Replace first with eighth
                meshmat_final{groups_final{i}(2), 5} = 'rx+rz';  % Replace second with third
                meshmat_final{groups_final{i}(3), 5} = 'rx+rz';  % Replace fourth with fifth
                meshmat_final{groups_final{i}(4), 5} = 'ry+rz'; % Replace sixth with seventh
            end
        end



        %% Append +z restraint to top nodes

        % check if all values of array are the same, in which case this
        % is a W_Y case, therefore the idx =2
        if numel(unique(cellfun(@(x) x, meshmat_final(:, 1)))) == 1
            idx=2;
        else
            idx=1;
        end

        init_restraint_sorted_data= sortrows(meshmat_final, [idx, 3]);

        % Find unique values in the idx column
        unique_first_column = unique(cell2mat(init_restraint_sorted_data(:, idx)));

        % Initialize the init_load_nod vector
        init_restraint_nod = {};

        % Iterate over unique values in the first column
        for i = 1:length(unique_first_column)
            % Find rows corresponding to the current unique value in the idx column
            rows = find(cell2mat(init_restraint_sorted_data(:, idx)) == unique_first_column(i));

            % Extract the subset of data for the current group (of same x node values)
            group_data = init_restraint_sorted_data(rows, :);

            % Find the maximum value(s) in the second column y for the current group
            max_values = max(cell2mat(group_data(:, 2)));

            % Find the rows where the second column matches the maximum value(s)
            max_rows = ismember(cell2mat(group_data(:, 2)), max_values);

            % Store the corresponding value from the sixth column in init_load_nod
            init_restraint_nod = [init_restraint_nod; group_data(max_rows, 6)];
        end


        for i = 1:length(meshmat_final(:, end))
            % Check if the current element of 7th col is in init_restraint_nod
            if any(strcmp(meshmat_final{i, 6}, init_restraint_nod))
                % If it is, append '+z' to the corresponding element in the 6th column
                meshmat_final{i, 5} = append(meshmat_final{i, 5}, '+z');
            end
        end

        %% Append x+y+z+rx+ry+rz restraint to bottom nodes

        % check if all values of array are the same, in which case this
        % is a W_Y case, therefore the idx =2
        if numel(unique(cellfun(@(x) x, meshmat_final(:, 1)))) == 1
            idx=2;
        else
            idx=1;
        end

        init_restraint_sorted_data= sortrows(meshmat_final, [idx, 3]);

        % Find unique values in the idx column
        unique_first_column = unique(cell2mat(init_restraint_sorted_data(:, idx)));

        % Initialize the init_restraint_nod vector
        init_restraint_nod = {};

        % Iterate over unique values in the first column
        for i = 1:length(unique_first_column)
            % Find rows corresponding to the current unique value in the idx column
            rows = find(cell2mat(init_restraint_sorted_data(:, idx)) == unique_first_column(i));

            % Extract the subset of data for the current group (of same x node values)
            group_data = init_restraint_sorted_data(rows, :);

            % Find the minimum value(s) in the second column y for the current group
            min_values = min(cell2mat(group_data(:, 2)));

            % Find the rows where the second column matches the minimum value(s)
            min_rows = ismember(cell2mat(group_data(:, 2)), min_values);

            % Store the corresponding value from the sixth column in init_restraint_nod
            init_restraint_nod = [init_restraint_nod; group_data(min_rows, 6)];
        end


        for i = 1:length(meshmat_final(:, end))
            % Check if the current element of 7th col is in init_restraint_nod
            if any(strcmp(meshmat_final{i, 6}, init_restraint_nod)) && ...
                    strcmp(meshmat_final{i, 5}, 'ry+rz')
                % If it is, append '+z' to the corresponding element in the 6th column
                meshmat_final{i, 5} = 'x+y+z+rx+ry+rz';
            end
        end

        %interchange the first and second row o bottom surface, to adapt to ADAPTIC syntax
        %ground level
        %         if strcmp(meshmat_final{1, 5}, 'x+y+z+rx+ry+rz')||...
        %                 strcmp(meshmat_final{1, 5}, 'ry+rz')
        %             meshmat_final{1, 5}='rx+rz';
        %             meshmat_final{2, 5}='x+y+z+rx+ry+rz';
        %         end



        %top level

        %         init_restraint_sorted_data= sortrows(meshmat_final, [idx, 3]);
        %
        %         % Find unique values in the idx column
        %         unique_first_column = unique(cell2mat(init_restraint_sorted_data(:, idx)));
        %
        %         % Initialize the init_restraint_nod vector
        %         init_restraint_nod = {};
        %
        %         % Iterate over unique values in the first column
        %         for i = 1:length(unique_first_column)
        %             % Find rows corresponding to the current unique value in the idx column
        %             rows = find(cell2mat(init_restraint_sorted_data(:, idx)) == unique_first_column(i));
        %
        %             % Extract the subset of data for the current group (of same x node values)
        %             group_data = init_restraint_sorted_data(rows, :);
        %
        %             % Find the minimum value(s) in the second column y for the current group
        %             max_values = max(cell2mat(group_data(:, 2)));
        %
        %             % Find the rows where the second column matches the minimum value(s)
        %             max_rows = ismember(cell2mat(group_data(:, 2)), max_values);
        %
        %             % Store the corresponding value from the sixth column in init_restraint_nod
        %             init_restraint_nod = [init_restraint_nod; group_data(max_rows, 6)];
        %         end
        %
        %         for i = 1:size(init_restraint_nod, 1)
        %
        %             if any(strcmp(meshmat_final{i, 6}, init_restraint_nod)) && ...
        %                     strcmp(meshmat_final{i, 5}, 'ry+rz')
        %                 % If it is, append '+z' to the corresponding element in the 6th column
        %                 meshmat_final{i, 5} = 'x+y+z+rx+ry+rz';
        %             end
        %         end


        %% REDEFINE CONNECTIVITY

        % Iterate through each 1x8 cell in elements_with_nodes
        for i = 1:size(elements_with_nodes, 1)
            % Get the element ID
            element_id = elements_with_nodes{i, 1};

            % Extract the 1x8 cell
            cell_data = elements_with_nodes{i, 3};

            % Get the first and last pairs of strings
            idx1 = find(strcmp(meshmat_final(:, 6), cell_data{1}));
            idx2 = find(strcmp(meshmat_final(:, 6), cell_data{end}));

            % Swap the strings if needed
            if strcmp(meshmat_final{idx1(1), 5}, 'rx+rz')||...
                    strcmp(meshmat_final{idx1(1), 5}, 'rx+rz+z')
                temp = cell_data{1};
                cell_data{1} = cell_data{end};
                cell_data{end} = temp;
            end

            loc=0;
            % Iterate through the rest of the pairs to reorder based on
            % restraints
            for j = 2:2:numel(cell_data)-2

                loc=loc+1; %every 2 nodes
                % Get the corresponding rows in meshmat_final
                idx1 = find(strcmp(meshmat_final(:, 6), cell_data{j}));
                idx2 = find(strcmp(meshmat_final(:, 6), cell_data{j+1}));

                % Swap the strings based on the order of 'ry+rz' and 'rx+rz'
                if loc==1 && ...
                        (strcmp(meshmat_final{idx1(1), 5}, 'rx+rz')||...
                        strcmp(meshmat_final{idx1(1), 5}, 'rx+rz+z'))
                    temp = cell_data{j};
                    cell_data{j} = cell_data{j+1};
                    cell_data{j+1} = temp;

                elseif loc==2 && ...
                        (strcmp(meshmat_final{idx1(1), 5}, 'ry+rz')||...
                        strcmp(meshmat_final{idx1(1), 5}, 'ry+rz+z')||...
                        strcmp(meshmat_final{idx1(1), 5}, 'x+y+z+rx+ry+rz'))
                    temp = cell_data{j};
                    cell_data{j} = cell_data{j+1};
                    cell_data{j+1} = temp;
                elseif loc==3 && ...
                        (strcmp(meshmat_final{idx1(1), 5}, 'rx+rz')||...
                        strcmp(meshmat_final{idx1(1), 5}, 'rx+rz+z'))
                    temp = cell_data{j};
                    cell_data{j} = cell_data{j+1};
                    cell_data{j+1} = temp;
                end

                % Update the 1x8 cell in elements_with_nodes
                elements_with_nodes{i, 3} = cell_data;
            end


            %% rearrange connectivity syntax


            % Initialize a cell array to store node IDs for each row in elem_with_nodes
            collected_IDs_all = cell(size(elements_with_nodes, 1), 1);
            
            % Iterate through each row of elem_with_nodes
            for i = 1:size(elements_with_nodes, 1)
                % Extract the 1x8 cell containing node IDs
                nodes_1x8 = elements_with_nodes{i, 3};
            
                % Initialize a cell array to store node IDs for the current row
                collected_IDs = cell(1, numel(nodes_1x8));
            
                % Iterate through each node ID in the 1x8 cell
                for j = 1:numel(nodes_1x8)
                    % Extract node ID from nodes_1x8{j}
                    node_id = nodes_1x8{j};
            
                    % Find the corresponding row in meshmat_final
                    idx = find(strcmp(meshmat_final(:, 6), node_id));
            
                    % Extract x and y coordinates from meshmat_final
                    node_x = meshmat_final{idx, 1};
                    node_y = meshmat_final{idx, 2};
            
                    % Find corresponding node IDs in meshmat_final with the same coordinates
                    idx = find(strcmp(cellfun(@num2str, meshmat_final(:, 1), 'UniformOutput', false), num2str(node_x))...
                        & strcmp(cellfun(@num2str, meshmat_final(:, 2), 'UniformOutput', false), num2str(node_y)));
            
                    % Store these corresponding node IDs in collected_IDs
                    collected_IDs{j} = meshmat_final(idx, 6);
                end
            
                % Store collected_IDs for the current row in collected_IDs_all
                collected_IDs_all{i} = collected_IDs;
            end
            
            % Initialize a cell array to store rearranged node IDs for each row in elem_with_nodes
            rearranged_nodeIDs = cell(size(elements_with_nodes, 1), 1);
            
            % Iterate through each row of elem_with_nodes
            for i = 1:size(elements_with_nodes, 1)
                % Extract the 1x8 cell containing collected node IDs
                collected_IDs = collected_IDs_all{i};
            
                % Initialize a cell array to store rearranged node IDs for the current row
                rearranged_nodeIDs_row = elements_with_nodes{i, 3};
                %cell(1, numel(collected_IDs));
            
                % Iterate through each node in the 1x8 cell
%                 for j = 1:numel(collected_IDs)
%                     % Assign node ID based on the specified conditions
%                     if j == 1 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'x+y+z+rx+ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'ry+rz+z'))
%             
%                         rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
% 
% 
%                     elseif j == 2 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'x+y+z+rx+ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'ry+rz+z'))
%                         rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
%                     elseif j == 3
%                         if numel(collected_IDs{j}) == 2 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'rx+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'rx+rz+z'))
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{end};
%                         elseif numel(collected_IDs{j}) == 3
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
%                         elseif numel(collected_IDs{j}) == 4
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{3};
%                         end
%                     elseif j == 4
%                         if numel(collected_IDs{j}) == 2 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'rx+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'rx+rz+z'))
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
%                         elseif numel(collected_IDs{j}) == 3
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
%                         elseif numel(collected_IDs{j}) == 4
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
%                         end
%                     elseif j == 5 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'x+y+z+rx+ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'ry+rz+z'))
%                         rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
%                     elseif j == 6 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'x+y+z+rx+ry+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{end})), 5}, 'ry+rz+z'))
%                         rearranged_nodeIDs_row{j} = collected_IDs{j}{end};
%                     elseif j == 7
%                         if numel(collected_IDs{j}) == 2 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'rx+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'rx+rz+z'))
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
%                         elseif numel(collected_IDs{j}) == 3
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
%                         elseif numel(collected_IDs{j}) == 4
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
%                         end
%                     elseif j==8
%                         if numel(collected_IDs{j}) == 2 && ...
%                             (strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'rx+rz')||...
%                             strcmp(meshmat_final{find(strcmp(meshmat_final(:, 6), collected_IDs{j}{1})), 5}, 'rx+rz+z'))
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
%                         elseif numel(collected_IDs{j}) == 3
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
%                         elseif numel(collected_IDs{j}) == 4
%                             rearranged_nodeIDs_row{j} = collected_IDs{j}{3};
%                         end
%                     end
%                 end
            
 % Iterate through each node in the 1x8 cell
                for j = 1:numel(collected_IDs)
                    % Assign node ID based on the specified conditions
                    if j == 1
                        rearranged_nodeIDs_row{j} = collected_IDs{j}{end};
                    elseif j == 2 
                        rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
                    elseif j == 3
                        if numel(collected_IDs{j}) == 2
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{end};
                        elseif numel(collected_IDs{j}) == 3
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{3};
                        elseif numel(collected_IDs{j}) == 4
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{3};
                        end
                    elseif j == 4
                        if numel(collected_IDs{j}) == 2 
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
                        elseif numel(collected_IDs{j}) == 3
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
                        elseif numel(collected_IDs{j}) == 4
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
                        end
                    elseif j == 5 
                        rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
                    elseif j == 6
                        rearranged_nodeIDs_row{j} = collected_IDs{j}{end};
                    elseif j == 7
                        if numel(collected_IDs{j}) == 2 
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
                        elseif numel(collected_IDs{j}) == 3
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
                        elseif numel(collected_IDs{j}) == 4
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
                        end
                    elseif j==8
                        if numel(collected_IDs{j}) == 2 
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{1};
                        elseif numel(collected_IDs{j}) == 3
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{2};
                        elseif numel(collected_IDs{j}) == 4
                            rearranged_nodeIDs_row{j} = collected_IDs{j}{3};
                        end
                    end
                end

                % Store rearranged node IDs for the current row in rearranged_nodeIDs
                rearranged_nodeIDs{i} = rearranged_nodeIDs_row;
            end
            
            elements_with_nodes(:,3)=rearranged_nodeIDs;


            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% INITIAL LOADS %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % check if all values of array are the same, in which case this
            % is a W_Y case, therefore the idx =2
            if numel(unique(cellfun(@(x) x, meshmat_final(:, 1)))) == 1
                idx=2;
            else
                idx=1;
            end

            init_load_sorted_data= sortrows(meshmat_final, [idx, 3]);

            % Find unique values in the idx column
            unique_first_column = unique(cell2mat(init_load_sorted_data(:, idx)));

            % Initialize the init_load_nod vector
            init_load_nod = {};
            %cell(length(unique_first_column), 1);

            % Iterate over unique values in the first column
            for i = 1:length(unique_first_column)
                % Find rows corresponding to the current unique value in the idx column
                rows = find(cell2mat(init_load_sorted_data(:, idx)) == unique_first_column(i));

                % Extract the subset of data for the current group (of same x node values)
                group_data = init_load_sorted_data(rows, :);

                % Find the maximum value in the second column y for the current group
                %[~, index] = max(cell2mat(group_data(:, 2)));

                % Store the corresponding value from the sixth column in init_load_nod
                %init_load_nod{i} = group_data{index, 6};

                % Find the maximum value(s) in the second column y for the current group
                max_values = max(cell2mat(group_data(:, 2)));

                % Find the rows where the second column matches the minimum value(s)
                max_rows = ismember(cell2mat(group_data(:, 2)), max_values);

                % Filter out elements based on corresponding strings in meshmat_final
                filter_indices = cellfun(@(x) strcmp(x, 'ry+rz') ||...
                    strcmp(x, 'ry+rz+z') ||...
                    strcmp(x, 'x+y+z+rx+ry+rz'), group_data(:, 5));

                % Store the corresponding value from the sixth column in init_load_nod
                init_load_nod = [init_load_nod; group_data(max_rows & filter_indices, 6)];

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


            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% MASS ELEMENTS %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            mass_elems_nodes=init_load_nod;
            mass_elems={};

            % Create a cell array of mass elements
            for i = 1:numel(mass_elems_nodes)
                mass_elems{i} = ['m' num2str(i)];
            end

            % Display the result
            mass_elems=mass_elems.';

            mass_elems_idx=repmat({'m1'}, 1, length(mass_elems_nodes)).';

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% TOP NODES VERT DISPLACEMENT %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            vert_dis_idx=size(mass_elems, 1)-1;

            % Initialize an empty cell array to store the elements
            vert_nodes = {};

            % Initialize the vector 'cpy'
            cpy_vec = cell(length(vert_dis_idx), 1).';

            % Generate elements starting with '10' and ending with digits 1 to 3
            for i = 1:vert_dis_idx
                % Create the element by concatenating '10' with the current digit
                node = ['10' num2str(i)];

                % Append the element to the cell array
                vert_nodes = [vert_nodes; node];

                % Create the element by repeating the string 'A' 'i' times
                cpy_vec{i,1} = 'cpy';

            end

            vert_dis_1=mass_elems_nodes(1:end-1);
            vert_dis_2=mass_elems_nodes(2:end);

            % g_forc_x vector

            gforcvec=mass_elems_nodes';


            %% SPRD information

            sprd_nodes= numel(gforcvec);
            sprd_ratios=[sprd_nodes, ones(1, sprd_nodes)];

            %% elem initial loads
            % Initialize a cell array to store element IDs
            in_lod_els = {};

            % Iterate through each row in elem_with_nodes
            for i = 1:size(elements_with_nodes, 1)
                % Check if the second column contains only 'macro1'
                if all(strcmp(elements_with_nodes{i, 2}, 'macro1'))
                    % Store the corresponding element ID
                    in_lod_els = [in_lod_els; elements_with_nodes{i, 1}];
                end
            end

            in_lod_typ=repmat({'udl1'}, 1, length(in_lod_els)).';
            in_lod_val=repmat({'0 -1.8e-5 0'}, 1, length(in_lod_els)).';






            %%%%%%%%%%%%%%%%%%%%
            %%%% SAVE FILES %%%%
            %%%%%%%%%%%%%%%%%%%%
            save_file_disc_corr(prefix, save_file_path, X_final*1000, Y_final*1000, Z_final*1000, node_names_final, 'sorted_coordinates');
            save_file_disc_corr(prefix, save_file_path, num_elements, elements_with_nodes, 'elem_coordinates');
            save_file_disc_corr(prefix, save_file_path, meshmat_final, 'restraints');
            save_file_disc_corr(prefix, save_file_path, init_load_nod, direction, type, value, 'initial_loads');
            save_file_disc_corr(prefix, save_file_path, mass_elems, mass_elems_idx, mass_elems_nodes, 'mass_elems');
            save_file_disc_corr(prefix, save_file_path, vert_nodes, cpy_vec, vert_dis_1, vert_dis_2, 'vert_disp');
            save_file_disc_corr(prefix, save_file_path, gforcvec, 'gforcvec');
            save_file_disc_corr(prefix, save_file_path, sprd_ratios, 'sprd');
            save_file_disc_corr(prefix, save_file_path, in_lod_els, in_lod_typ, in_lod_val, 'in_lod');



        end
        %% Numbers to be inputted in plot_deformed_shape

        disp([prefix, 'nbn = ', num2str(size(meshmat_final, 1))]);
        disp([prefix, 'nbe = ', num2str(size(elements_with_nodes, 1))]);
    end
    successOutput = true;

end
