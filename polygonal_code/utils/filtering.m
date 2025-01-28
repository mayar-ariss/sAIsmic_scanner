
% Mayar Ariss 25 Mar 2024

clear; clc; close all;

%%% SPANDRELS %%%

% Preallocate surface areas array
num_surfaces_spandrels = size(XSpandrels, 2);
surfaceAreas_spandrels = zeros(1, num_surfaces_spandrels);

spandrel_diag=zeros(size(surfaceAreas_spandrels));
vertices1=zeros(size(surfaceAreas_spandrels));
vertices2=zeros(size(surfaceAreas_spandrels));

XSpandrels_filtered=[];
YSpandrels_filtered=[];
ZSpandrels_filtered=[];
j_spandrel=0;

% Calculate surface areas for each surface
for i = 1:num_surfaces_spandrels

    spandrel_corner1=[XSpandrels(1,i), YSpandrels(1,i), ZSpandrels(1,i)];
    spandrel_corner2=[XSpandrels(2,i), YSpandrels(2,i), ZSpandrels(2,i)];
    spandrel_corner3=[XSpandrels(3,i), YSpandrels(3,i), ZSpandrels(3,i)];
    spandrel_corner4=[XSpandrels(4,i), YSpandrels(4,i), ZSpandrels(4,i)];

    spandrel_diag(i) = sqrt(sum((permute(spandrel_corner1, [1 3 2]) - permute(spandrel_corner3, [3 1 2])) .^ 2, 3));
    spandrel_vertices1(i)=sqrt(sum((permute(spandrel_corner1, [1 3 2]) - permute(spandrel_corner2, [3 1 2])) .^ 2, 3));
    spandrel_vertices2(i)=sqrt(sum((permute(spandrel_corner1, [1 3 2]) - permute(spandrel_corner4, [3 1 2])) .^ 2, 3));

    if spandrel_diag(i) == spandrel_vertices1(i) || spandrel_diag(i) == spandrel_vertices2(i)

    else
        j_spandrel=j_spandrel+1;
        XSpandrels_filtered(:,j_spandrel)=XSpandrels(:,i);
        YSpandrels_filtered(:,j_spandrel)=YSpandrels(:,i);
        ZSpandrels_filtered(:,j_spandrel)=ZSpandrels(:,i);

    end
    
end

%%% PIERS %%%

num_surfaces_piers = size(XPiers, 2);
surfaceAreas_piers = zeros(1, num_surfaces_piers);

pier_diag=zeros(size(surfaceAreas_piers));

XPiers_filtered=[];
YPiers_filtered=[];
ZPiers_filtered=[];
j_pier=0;

% Calculate surface areas for each surface
for i = 1:num_surfaces_piers

    pier_corner1=[XPiers(1,i), YPiers(1,i), ZPiers(1,i)];
    pier_corner2=[XPiers(2,i), YPiers(2,i), ZPiers(2,i)];
    pier_corner3=[XPiers(3,i), YPiers(3,i), ZPiers(3,i)];
    pier_corner4=[XPiers(4,i), YPiers(4,i), ZPiers(4,i)];

    pier_diag(i) = sqrt(sum((permute(pier_corner1, [1 3 2]) - permute(pier_corner3, [3 1 2])) .^ 2, 3));
    pier_vertices1(i)=sqrt(sum((permute(pier_corner1, [1 3 2]) - permute(pier_corner2, [3 1 2])) .^ 2, 3));
    pier_vertices2(i)=sqrt(sum((permute(pier_corner1, [1 3 2]) - permute(pier_corner4, [3 1 2])) .^ 2, 3));

    if pier_diag(i) == pier_vertices1(i) || pier_diag(i) == pier_vertices2(i)

    else
        j_pier=j_pier+1;
        XPiers_filtered(:,j_pier)=XPiers(:,i);
        YPiers_filtered(:,j_pier)=YPiers(:,i);
        ZPiers_filtered(:,j_pier)=ZPiers(:,i);

    end
    
end

%%% NODES %%%

num_surfaces_nodes = size(XNodes, 2);
surfaceAreas_nodes = zeros(1, num_surfaces_nodes);

node_diag=zeros(size(surfaceAreas_nodes));

XNodes_filtered=[];
YNodes_filtered=[];
ZNodes_filtered=[];
j_node=0;

% Calculate surface areas for each surface
for i = 1:num_surfaces_nodes

    node_corner1=[XNodes(1,i), YNodes(1,i), ZNodes(1,i)];
    node_corner2=[XNodes(2,i), YNodes(2,i), ZNodes(2,i)];
    node_corner3=[XNodes(3,i), YNodes(3,i), ZNodes(3,i)];
    node_corner4=[XNodes(4,i), YNodes(4,i), ZNodes(4,i)];

    node_diag(i) = sqrt(sum((permute(node_corner1, [1 3 2]) - permute(node_corner3, [3 1 2])) .^ 2, 3));
    node_vertices1(i)=sqrt(sum((permute(node_corner1, [1 3 2]) - permute(node_corner2, [3 1 2])) .^ 2, 3));
    node_vertices2(i)=sqrt(sum((permute(node_corner1, [1 3 2]) - permute(node_corner4, [3 1 2])) .^ 2, 3));

    if node_diag(i) == node_vertices1(i) || node_diag(i) == node_vertices2(i)

    else
        j_node=j_node+1;
        XNodes_filtered(:,j_node)=XNodes(:,i);
        YNodes_filtered(:,j_node)=YNodes(:,i);
        ZNodes_filtered(:,j_node)=ZNodes(:,i);

    end
    
end
