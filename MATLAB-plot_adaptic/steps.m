function [ nb_steps ] = steps( file,name )
% Gets the total number of steps stored in file struc#file.num

if nargin < 2
    name='struc';
end

if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

nb_steps=0;
oneline=fgets(fid);
while ischar(oneline)
    oneline=fgets(fid);
  
    % Step count
    if length(oneline)>3 && strcmp(oneline(1:4),'#io1')
        nb_steps=nb_steps+1
    end
end
fclose(fid);
end

