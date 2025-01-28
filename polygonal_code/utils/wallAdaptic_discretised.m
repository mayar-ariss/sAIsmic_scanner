% Mayar Ariss 10 May 2024

function [successOutput, node_name] = wallAdaptic_discretised(XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face)


% Combine X_face and Y_face into a single cell array
faces = {X_face, Y_face};

for face_idx = 1:length(faces)
    face_values = faces{face_idx};
    for i = 1:length(face_values)
        values = face_values(i);
        zero_indices = [];


        %% NODES
        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XNodes_filtered, 3) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YNodes_filtered, 3) == values);
        end

        % Select X, Y, and Z for only front frame
        X_nodes = round(XNodes_filtered(:, zero_indices), 3);
        Y_nodes = round(YNodes_filtered(:, zero_indices), 3);
        Z_nodes = round(ZNodes_filtered(:, zero_indices), 3);

        %% SPANDRELS

        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XSpandrels_filtered, 3) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YSpandrels_filtered, 3) == values);
        end

        % Select X, Y, and Z for only front frame
        X_spandrels = round(XSpandrels_filtered(:, zero_indices), 3);
        Y_spandrels = round(YSpandrels_filtered(:, zero_indices), 3);
        Z_spandrels = round(ZSpandrels_filtered(:, zero_indices), 3);


        %% PIERS

        if face_idx == 1
            % X_face processing
            zero_indices = all(round(XPiers_filtered, 3) == values);
        elseif face_idx == 2
            % Y_face processing
            zero_indices = all(round(YPiers_filtered, 3) == values);
        end

        X_piers = round(XPiers_filtered(:, zero_indices), 3);
        Y_piers = round(YPiers_filtered(:, zero_indices), 3);
        Z_piers = round(ZPiers_filtered(:, zero_indices), 3);



        %% BLOCKS

        % NODE BLOCKS
        % Number of NODE block elements
        num_node_blocks = size(X_nodes, 2); % Assuming X_nodes and Y_nodes have the same number of columns       
        % Initialize block variable to store coordinates
        node_blocks = zeros(4, 2, num_node_blocks); % 4 rows for x and y coordinates, 2 columns for x and y
       
        % Iterate over each node block element
        for i = 1:num_node_blocks
            % Extract x and y coordinates for the current block element
            x_coords = X_nodes(:, i);
            y_coords = Y_nodes(:, i);
            z_coords = Z_nodes(:, i);

            % Store coordinates in the blocks variable
            if numel(unique(X_nodes(:))) == 1
                node_blocks(:,:,i) = [y_coords, z_coords];
            elseif numel(unique(Y_nodes(:))) == 1
                node_blocks(:,:,i) = [x_coords, z_coords];
            end

        end


        % SPANDREL BLOCKS
        % Number of SPANDREL block elements
        num_spdrl_blocks = size(X_spandrels, 2); % Assuming X_spandrels and Y_spandrels have the same number of columns       
        % Initialize block variable to store coordinates
        spdrl_blocks = zeros(4, 2, num_spdrl_blocks); % 4 rows for x and y coordinates, 2 columns for x and y
       
        % Iterate over each spandrel block element
        for i = 1:num_spdrl_blocks
            % Extract x and y coordinates for the current block element
            x_coords = X_spandrels(:, i);
            y_coords = Y_spandrels(:, i);
            z_coords = Z_spandrels(:, i);

            % Store coordinates in the blocks variable
            if numel(unique(X_spandrels(:))) == 1
                spdrl_blocks(:,:,i) = [y_coords, z_coords];
            elseif numel(unique(Y_spandrels(:))) == 1
                spdrl_blocks(:,:,i) = [x_coords, z_coords];
            end

        end

        % PIER BLOCKS
        % Number of PIER block elements
        num_pier_blocks = size(X_piers, 2); % Assuming X_piers and Y_piers have the same number of columns       
        % Initialize block variable to store coordinates
        pier_blocks = zeros(4, 2, num_pier_blocks); % 4 rows for x and y coordinates, 2 columns for x and y
       
        % Iterate over each pier block element
        for i = 1:num_pier_blocks
            % Extract x and y coordinates for the current block element
            x_coords = X_piers(:, i);
            y_coords = Y_piers(:, i);
            z_coords = Z_piers(:, i);

            % Store coordinates in the blocks variable
            if numel(unique(X_piers(:))) == 1
                pier_blocks(:,:,i) = [y_coords, z_coords];
            elseif numel(unique(Y_piers(:))) == 1
                pier_blocks(:,:,i) = [x_coords, z_coords];
            end

        end



        %% DISCRETISATION

        % Initialize new_blocks with the original blocks
        new_blocks = pier_blocks;

        for i=1:num_pier_blocks

            % Select the coordinates of the first block
            my_block = pier_blocks(:,:,i);

            % Define the coordinates of all other blocks
            other_blocks = pier_blocks;
            % Remove the current block
            other_blocks(:,:,i) = [];

            % Number of block elements
            num_other_blocks = size(other_blocks, 3); % Assuming X_nodes and Y_nodes have the same number of columns

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
            for j = 1:num_other_blocks
                
                other_block = other_blocks(:,:,j);

                % Check if any vertex of any other block lies on any side of the current block
                for k = 1:4
                    vertex = other_block(k,:);
                    for side = 1:4
                        my_side = my_block([side, mod(side, 4) + 1], :);

                        if vertex(1,1)==my_side(1,1) && vertex(1,1)== my_side(2,1) ...
                                && min(my_side(:,2)) < vertex(1,2) && vertex(1,2) < max(my_side(:,2))
                            disp(['FACE ', num2str(values),' Vertex ', num2str(k), ' of my block ', num2str(i), ' lies on side ', num2str(side), ' of other block ', num2str(j)]);
                        new_block1= [sort(my_block(:, 1)), [min(my_block(:,1)), vertex(1,2), vertex(1,2), min(my_block(:,1))]'];
                        new_block2= [sort(my_block(:, 1)), [vertex(1,2), max(my_block(:,1)), max(my_block(:,1)), vertex(1,2)]'];
                        
                        % Append new blocks to new_blocks
                        new_blocks = cat(3, new_blocks, new_block1, new_block2);    

                        % Remove my_block from new_blocks
                        new_blocks(:,:,i) = [];

                        end

                    end
                end

            end


        end

    end

end

end
