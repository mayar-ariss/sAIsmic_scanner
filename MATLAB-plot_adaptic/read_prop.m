function [ SPRP ] = read_prop( name,file )
% Gets the nodal names (NN) and coordinates (XYZ) from the input .dat file

% nbn = number of nodes to be read
% name = name of the input file name.dat or name#00i.dat
% file = number of partition (0 for parent file)

if nargin < 2
    file=0;
end

SPRP=zeros(5,1);

if file>0
    fid=fopen(sprintf('../%01s#%03d.dat',name,file),'r');
else
    fid=fopen(sprintf('../%01s.dat',name),'r');
end

oneline=fgets(fid);
inod=0;
fin=false;

while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the structural.nodal 
    if ~isempty(strfind(oneline,'matd')) && ~isempty(strfind(oneline,'mcrd'))
        oneline(strfind(oneline,'matd'):strfind(oneline,'matd')+3)=[];
        oneline(strfind(oneline,'mcrd'):strfind(oneline,'mcrd')+3)=[];
        oneline(strfind(oneline,'&'))=[];
        node=str2num(oneline);  SPRP(1:4)=node;
        %skip lines containing headers or comments
		oneline=fgets(fid);
		while ~isempty(strfind(oneline,'#'))
			oneline=fgets(fid);
		end
		oneline(strfind(oneline,'&'))=[];
        node=str2num(oneline);  SPRP(5:10)=node;
		oneline=fgets(fid);
		while ~isempty(strfind(oneline,'#'))
			oneline=fgets(fid);
		end
		oneline(strfind(oneline,'&'))=[];
        node=str2num(oneline);  SPRP(11:19)=node;
		oneline=fgets(fid);
		while ~isempty(strfind(oneline,'#'))
			oneline=fgets(fid);
		end
		node=str2num(oneline);  SPRP(20:23)=node;
        break
    end
end
fclose(fid);
end
