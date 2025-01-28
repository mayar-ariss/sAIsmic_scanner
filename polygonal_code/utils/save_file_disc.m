% Mayar Ariss 06/05/2024

function save_file_disc(prefix, base_path, varargin)
    output_file = base_path;
    [path, name, ext] = fileparts(base_path);
    output_file = fullfile(path, prefix + "_" + varargin{end} + ext);

    % Combine all arguments except the last one
    data = varargin(1:end-1);
    % Extract the identifier string
    identifier = varargin{end};
    % Open the file for writing
    fid = fopen(output_file, 'w');
    if fid == -1
        error('Could not open file for writing: %s', output_path);
    end
    
    % Write data to the file based on the identifier
    switch identifier
        case 'sorted_coordinates'
            % Write header
            fprintf(fid, 'nod.name\tx\ty\tz\n');
            % Write data
            for i = 1:numel(data{1})
                fprintf(fid, '%s\t%f\t%f\t%f\n', data{4}{i}, data{1}(i), data{2}(i), data{3}(i));
            end
        case 'elem_coordinates'
            % Write header
            fprintf(fid, 'elm.name\tgrp.name \tnod.name\n');
            % Write data
            for i = 1:data{1}
                fprintf(fid, '%s\t%s\t%s\n', data{2}{i, 1}, data{2}{i, 2}, strjoin(data{2}{i, 3}, ' '));
            end
        case 'restraints'
            % Write header
            fprintf(fid, 'nod.name\tdirection\n');
            % Write data
            for i = 1:size(data{1}, 1)
                fprintf(fid, '%s\t%s\n', data{1}{i, 6}, data{1}{i, 5}{1});
            end
        case 'initial_loads'
            % Write header
            fprintf(fid, 'nod.name\tdirection\ttype\tvalue\n');
            % Write data
            for i = 1:numel(data{1})
                fprintf(fid, '%s\t%s\t%s\t%d\n', data{1}{i}, data{2}{i}, data{3}{i}, data{4}(i));
            end
        otherwise
            error('Invalid identifier: %s', identifier);
    end
    
    % Close the file
    fclose(fid);

end
