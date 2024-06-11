function [ ] = plot12dofBeam( pos,rot,beam,T3D1,T3D2,color )
% Plots a beam in 3D space considering position and rotations at edges.
% 3 corrdinates and 2 rotations per node are considered

% ATTENTION: length of beam is calculated based on current position of
% edges => different from small displacement assumption, but practically no
% difference

% Input:
% pos(i,1:3) = position of node i (i=1,2)
% rot(i,1:2) = roations of node i (i=1,2)

if nargin<5
    color=[0,0,0];
end

% transform coordinates to local system of the beam
if beam==1 || beam==3
    lpos=(T3D1*pos')';
    lrot=(T3D1*rot')';
else
    lpos=(T3D2*pos')';
    lrot=(T3D2*rot')';
end
dofs=[lpos(1,:),lrot(1,[2,3]),lpos(2,:),lrot(2,[2,3])]';

dist=lpos(1,1:3)-lpos(2,1:3);
L=sqrt(dot(dist,dist));

% local coordinate
xi=-1:0.2:1;

% position of each control point
lxyz=zeros(length(xi),3);
for p=1:length(xi)
    N=shapefunct_beam(xi(p),L);
    lxyz(p,:)=(N*dofs)';
end

% Transform back to global axis
if beam==1 || beam==3
    xyz=(T3D1'*lxyz')';
else
    xyz=(T3D2'*lxyz')';
end

plot3(xyz(:,1),xyz(:,2),xyz(:,3),'Linestyle','-','Color',color,'LineWidth',2)

end

