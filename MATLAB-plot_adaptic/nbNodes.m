function [ nb_nodes ] = nbNodes( file,name )
% Gets the total number of steps stored in file struc#file.num

if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

nb_nodes=0;
oneline=fgets(fid);
while ischar(oneline)
    oneline=fgets(fid);
  
    % Step count
    if ~isempty(strfind(oneline,'#in2'))
        oneline=fgets(fid);
        oneline=fgets(fid);
        oneline=fgets(fid); % first node
        while ~isempty(str2num(oneline))
            nb_nodes=nb_nodes+1;
            oneline=fgets(fid);
        end
        break
    end
end
fclose(fid);
end

