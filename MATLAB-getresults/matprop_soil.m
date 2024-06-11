function [ SOIL ] = matprop_soil( file )
% Gets the displacements of the input "node" for each step, from the output 
% file struc#file.num

% The output disp is a matrix and contains the displacements in x,y and z 
% of the node at each step

if file>0
    fid=fopen(sprintf('../struc#%03d.dat',file),'r');
else
    fid=fopen('../struc.dat','r');
end

oneline=fgets(fid);

while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the interface material properties
    if length(oneline)>4 && strcmp(oneline(1:7),' sblock')
        oneline=oneline(20:end-2);
        sblock=str2num(oneline)';
        oneline=fgets(fid);
        oneline=oneline(20:end-2);
        sblock=[sblock;str2num(oneline)'];
        oneline=fgets(fid);
        oneline=oneline(20:end);
        sblock=[sblock;str2num(oneline)'];
    end
    if length(oneline)>4 && strcmp(oneline(1:7),' swalls')
        oneline=oneline(20:end);
        swalls=str2num(oneline)';
        break
    end
    
end
fclose(fid);

SOIL=cell(size(sblock,1)+1,3);
SOIL{2,1}='Esoil';
SOIL{3,1}='v';
SOIL{4,1}='COH_in';
SOIL{5,1}='tanf_in';
SOIL{6,1}='COH_res';
SOIL{7,1}='tanf_res';
SOIL{9,1}='hyp_approx_tension';

SOIL{1,2}='sblock';
SOIL{1,3}='swalls';
for i=1:size(sblock,1)
    SOIL{i+1,2}=sblock(i);
end
for i=1:size(swalls,1)
    SOIL{i+1,3}=swalls(i);
end

end
