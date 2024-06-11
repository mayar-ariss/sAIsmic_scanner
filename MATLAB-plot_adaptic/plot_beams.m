function [ ] = plot_beams( NXYZ,rot,T3D1,T3D2,color )
% Plots the external beams of one macroelement

% Input:
% - NXYZ(8,1+ndof) = node numbers, current 3d positions and rotations of the 
%                    8 nodes of the macroelement in connectivity order
%                    ndof = 3 when plotting undeformed shape (we only need the positions)
%                         = 6 otherwise (positions and rotations), though some 
%                           of the rotations may not be active
% - rot = false(default)/true defines whether there are rotational dofs on the edges

if nargin<4
    color=[0,0,0];
end

if rot && size(NXYZ,2)<7
    sprintf(' Error: Rotations for the plotting of the beams are not provided')
    return
end
    
if ~rot
    pos=NXYZ(:,2:4);
    for beam=1:4
        plot3(pos([2*beam,2*beam-1],1),pos([2*beam,2*beam-1],2),pos([2*beam,2*beam-1],3),...
            'Linestyle','-','LineWidth',2,'Color',color)
    end
 else
     pos=NXYZ(:,2:4);
     rot=NXYZ(:,5:7);
     for beam=1:2
         plot12dofBeam(pos([2*beam-1,2*beam],:),rot([2*beam-1,2*beam],:),beam,T3D1,T3D2,color)
     end
     for beam=3:4
         plot12dofBeam(pos([2*beam,2*beam-1],:),rot([2*beam,2*beam-1],:),beam,T3D1,T3D2,color)
     end
 end

end

