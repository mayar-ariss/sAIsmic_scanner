function [] = dam_undeformed( NN,XYZ,CON,EN,DEN,DAM,thick,Vert)
% Plots the undeformed shape of the URM macroelements

% Input:
% - NN{pp}(nnod,1) = nodal names of all nodes at each partition pp
% - XYZ{pp}(nnod,3) = initial positions of all nodes at each partition pp
% - CON{pp}(nelem,8) = macroelement connectivity for each partition pp
% - DEN{pp}(ndamelem,1) = list of element names that their damage level will
%                          be plotted
% - DAMS{pp}(ndamelem,1) = damage of shear springs
% - DAMO{pp}(ndamelem,1) = damage of out-of-plane springs
% - thick = thickness of plotted macrelement blocks

beamcolor=[0,0,1];
blockcolor=[0,0,0];       

figure
hold on

% for each partition
for pp=1:length(CON)    
    if strcmp(DEN{pp},'all')
        DINCL=true(length(EN{pp}),1);
    else
        DINCL=ismember(EN{pp},DEN{pp});
    end

    % plot each macroelement's beams and inner blocks
    for ee=1:size(CON{pp},1)
        
        % obtain matrix with nodal names and initial positions of
        % macroelement nodes
        nxyz0=matrix_part(CON{pp}(ee,:),1,[NN{pp},XYZ{pp}],false);
        
        % undeformed shape -> no rotations (straight beams)
        rot=false;         
        
        % current adofs are all zero
        cadf=zeros(8,1);
        
        % Obtain transofrmation matrices (not used for undeformed)
        [~,T3D1,T3D2]=get_plane(nxyz0([1,2,4,6],2:end));
        
        plot_beams(nxyz0,rot,T3D1,T3D2,beamcolor);
        plot_block(nxyz0,cadf,thick,blockcolor,Vert);
        if DINCL(ee)
            plot_damage(nxyz0,cadf,DAM{pp}(ee));
        end
    end
    
end

end
