function [ CADOFS ] = get_CAdofs( file,name,nbe,nst,plotscale )
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

CADOFS=cell(length(file),1);

for pp=1:length(file)
    CADOFS{pp}=read_cadofs(file(pp),name,nbe(pp),nst);
    CADOFS{pp}=plotscale*CADOFS{pp};
end

end

