function [ LF ] = get_LF( file,name )
%Gets the load factor at each step from file struc#file.num

if nargin < 2
    name='struc';
end

if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

oneline=fgets(fid);
ilf=0;
while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the Load Factor (LF) of the current step (ilf)
    if length(oneline)>3 && strcmp(oneline(1:4),'#io1')
        ilf=ilf+1;
        oneline=fgets(fid);
        step=fgets(fid);
        step=str2num(step);
        LF(ilf,1)=step(3);
    end
end
fclose(fid);
end

