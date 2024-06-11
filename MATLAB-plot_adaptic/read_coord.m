function [ NN,XYZ ] = read_coord( nbn,name,file )
% Gets the nodal names (NN) and coordinates (XYZ) from the input .dat file

% nbn = number of nodes to be read
% name = name of the input file name.dat or name#00i.dat
% file = number of partition (0 for parent file)

if nargin < 2
    file=0;
end

NN=zeros(nbn,1);
XYZ=zeros(nbn,3);

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
    if ~isempty(strfind(oneline,'structural'))
        %skip lines containing headers or comments
        while ~fin
            oneline=fgets(fid);
            node=str2num(oneline);
            while ~isempty(strfind(oneline,'#')) ||...
                    sum(node)==0
                oneline=fgets(fid);
                node=str2num(oneline);
            end        
            inod=inod+1;
            NN(inod)=node(1);
            XYZ(inod,:)=node(2:end);
            if inod==nbn
            %if sum(node)==0
                fin=true;
            end
        end
        break
    end
end
fclose(fid);
end
