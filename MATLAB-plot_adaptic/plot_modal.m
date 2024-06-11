function [] = plot_modal( NN,XYZ0,CXYZ,CON,rot,nst )
%Plots the undeformed shape of the URM macroelements

% Input:
% - NN{pp}(nnod,1) = nodal names of all nodes at each partition pp
% - XYZ0{pp}(nnod,3) = initial positions of all nodes at each partition pp
% - CXYZ{pp}(nnod,3,nstep) = current positions (for every step) of all nodes at each partition pp
% - CON{pp}(nelem,8) = macroelement connectivity for each partition pp
% - rot{pp} = flag (t/f) for each partition defining whether rotations of external beams are plotted
% - CADOF{pp}(8*nelm,2) = numbers and current values of all additional dofs of each partition
% - ACON{pp}(nelm,8) = adof connectivity for each element
% - thick = thickness of plotted macrelement blocks

%================
% CADOF,ACON: maybe 1 matrix will be enough, depending on how the adofs are
% written in the num file:
% - if they are written as the displacements, then I keep the two matrices 
% - if they are written as element variables, then I will store in one
%   matrix directly
%================

beamcolor=[0,0,0];
%blockcolor=[0,0,0];        

figure(1)
clf

% for each partition
for pp=1:length(CON)
    
    % plot each macroelement's beams and inner blocks
    for ee=1:size(CON{pp},1)
        
        % obtain matrix with nodal names, current positions and rotations of 
        % macroelement nodes (to be used for plotting external beams)
        nxyzr=matrix_part(CON{pp}(ee,:),1,CXYZ{pp}(:,:,nst),false);
        
        % obtain matrix with nodal names, and initial positions of
        % macroelement nodes (to be used for plotting block)
        nxyz0=matrix_part(CON{pp}(ee,:),1,[NN{pp},XYZ0{pp}],false);
        [~,T3D1,T3D2]=get_plane(nxyz0([1,2,4,6],2:end));
        
        % obtain matrix with current values of additional dofs of
        % macroelement (to be used for plotting block)
        
        %view(2)
        figure(1)
        hold on
        %axis off
        set(gcf,'color','w')
        %axis equal
%        plot_block(nxyz0,CADOF{pp}(ee,:,nst),thick,blockcolor);
        plot_beams(nxyzr,rot{pp},T3D1,T3D2,beamcolor);        
        
    end
end
%pause(0.5)

end