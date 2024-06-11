function [ disp ] = read_modal( file,name,nmod )
% Obtain displacements and rotations for all nodes of the system at each
% step
% Input:
% - file = number of partition (0 for parent)
% - name = name of dat file
% Output:
% - disp(nnod,1+ndof,nstep) = matrix containing nodal name, translations and 
%                           rotations for each node at each step.
%                           nnod = total number of nodes
%                           ndof = number of dofs per node
%                           nstep = number of steps performed in the
%                           analysis

%nst=steps(file,name);
nnod=nbNodes(file,name);

disp=zeros(nnod,7,nmod);

if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

oneline=fgets(fid);
istep=0;
while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the displacement of the nodes in 'node' (disp{node}) at the current step (ilf)
    if ~isempty(strfind(oneline,'#in2'))
        istep=istep+1
        oneline=fgets(fid); oneline=fgets(fid); oneline=fgets(fid);
        inod=0;
        % Read every line of nodal displacement until you find all the
        % nodes in 'node'
        while inod<nnod
            inod=inod+1;
            disp(inod,:,istep)=str2num(oneline);
            oneline=fgets(fid);
        end
        if istep==nmod
            fclose(fid);
            break
        end
    end
end

end

