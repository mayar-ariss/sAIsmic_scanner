function [ DispRot,nst ] = get_DispRot( file,name )
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

nst=steps(0,name);
DispRot=cell(length(file),1);

for pp=1:length(file)
    DispRot{pp}=read_ur(file(pp),name,nst);
end

end

