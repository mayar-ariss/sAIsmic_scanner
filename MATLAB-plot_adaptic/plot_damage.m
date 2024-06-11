function [ ] = plot_damage( NXYZ,ADOF,DAM )
% Plots inner block of one macroelement

% Input:
% - NXYZ(8,1+3)= node numbers and initial 3d positions of the 8 nodes of the
%                macroelement in connectivity order
% - cadf = values of additional dofs of the element at each step of the
%          analysis
% - thick = thickness of plotted block
% - color = color vector for inner block lines

ind=[1,2,4,6];
% extract local displacements for each edge of the block
dlocal=[ADOF(1),-ADOF(4),ADOF(5);...
      ADOF(1),ADOF(2),ADOF(6);...
      -ADOF(3),ADOF(2),ADOF(7);...
      -ADOF(3),-ADOF(4),ADOF(8)];

dxyz=transform3d(NXYZ,dlocal);
  
pos=NXYZ(ind,2:4)+dxyz;

con=[1,2,3,4];

color=DAM;
caxis([0 1])
colormap(jet)
patch(pos(con,1),pos(con,2),pos(con,3),color)
colorbar

view(2)
hold on


% % normals to the surface at edges
% normal=zeros(4,3);
% normal(1,:)=cross(pos(4,:)-pos(1,:),pos(2,:)-pos(1,:));
% normal(2,:)=cross(pos(1,:)-pos(2,:),pos(3,:)-pos(2,:));
% normal(3,:)=cross(pos(2,:)-pos(3,:),pos(4,:)-pos(3,:));
% normal(4,:)=cross(pos(3,:)-pos(4,:),pos(1,:)-pos(4,:));
% for i=1:4
%     normal(i,:)=normal(i,:)/norm(normal(i,:));
% end
% 
% % front and back face
% posf=zeros(4,3);
% posb=zeros(4,3);
% for i=1:4
%     posf(i,:)=pos(i,:)+normal(i,:)*thick/2;
%     posb(i,:)=pos(i,:)-normal(i,:)*thick/2;
% end
% 
% plot3(posf(con,1),posf(con,2),posf(con,3),...
%         'Linestyle','-','Color',color)
% plot3(posb(con,1),posb(con,2),posb(con,3),...
%         'Linestyle','-','Color',color)
% for i=1:4
%     plot3([posf(con(i),1);posb(con(i),1)],...
%           [posf(con(i),2);posb(con(i),2)],...
%           [posf(con(i),3);posb(con(i),3)],...
%           'Linestyle','-','Color',color)
% end

