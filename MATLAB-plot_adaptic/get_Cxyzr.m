function [ Cxyzr ] = get_Cxyzr( NN,XYZ,UR,scale )
% For each partition containing macroelements, the function provides the 
% current positions and rotations of all the nodes of the system/partition 
% at each step of the analysis

% Input:
% - NN{npart}(mnnod,1) = nodal names of the macroelement nodes in the partition
% - XYZ{npart}(mnnod,3) = initial coordinates of the macroelement nodes in the partition
% - UR{npart}(nnod,7,nstep) = total displacements and total rotations of each node at each step
% - Scale factor for displacements and rotations for plotting reasons

% Output:
% - CXYZ{npart}(nnod,6,nstep) = current positions and total rotations of each macroelement node at each step

if nargin<3
    scale=1;
end
nst=size(UR{1},3);
Cxyzr=cell(size(XYZ,1),1);

for pp=1:size(XYZ,1)
    for ss=1:nst
        Cxyzr{pp}(:,:,ss) = matrix_part( NN{pp},1,UR{pp}(:,:,ss),false);
        Cxyzr{pp}(:,2:end,ss) = scale*Cxyzr{pp}(:,2:end,ss);
        Cxyzr{pp}(:,2:4,ss) = Cxyzr{pp}(:,2:4,ss) + XYZ{pp};        
    end
end

end

