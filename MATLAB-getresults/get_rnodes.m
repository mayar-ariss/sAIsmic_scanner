function [ nodes ] = get_rnodes( file,name )
% Obtains the nodes to which dynamic loads are applied

if nargin < 2
    name='struc';
end

% open .dat file
if file>0
    fid=fopen(sprintf('../%01s#%03d.dat',name,file),'r');
else
    fid=fopen(sprintf('../%01s.dat',name),'r');
end

% initiate counter of nodes stored
in=0;

% read file
oneline=fgets(fid);
while ischar(oneline)
    oneline=fgets(fid);
    if length(oneline)>1
        if ~isempty(strfind(oneline,'dynamic.loads')) || ~isempty(strfind(oneline,'proportional.loads'))
            oneline=fgets(fid);
            rnod=fgets(fid);
            if isempty(strfind(rnod,' x ')) && isempty(strfind(rnod,' y ')) && isempty(strfind(rnod,' z '))
                rnod=fgets(fid);
            end
            while ~isempty(strfind(rnod,'x')) || ~isempty(strfind(rnod,'y')) || ~isempty(strfind(rnod,'z'))
                in=in+1;
                % - find position of 'x' or 'y' or 'z', which is the
                %   direction input right after the node number
                % - keep only the part of the line before the direction
                %   input - this part contains the node number
                % - convert to number and store in 'nodes'
                if ~isempty(strfind(rnod,'x'))
                    rnod=rnod(1:strfind(rnod,'x')-1);
                    nodes(in)=str2double(rnod);                
                elseif ~isempty(strfind(rnod,'y'))
                    rnod=rnod(1:strfind(rnod,'y')-1);
                    nodes(in)=str2double(rnod);                   
                elseif ~isempty(strfind(rnod,'z'))
                    rnod=rnod(1:strfind(rnod,'z')-1);
                    nodes(in)=str2double(rnod);
                end
                rnod=fgets(fid);
            end
        end
    end
end
fclose(fid);
end

