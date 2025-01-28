


% Preallocate surface areas array
num_surfaces_spandrels = size(XSpandrels, 2);
surfaceAreas_spandrels = zeros(1, num_surfaces_spandrels);

diag=zeros(size(surfaceAreas_spandrels));
vertices1=zeros(size(surfaceAreas_spandrels));
vertices2=zeros(size(surfaceAreas_spandrels));

XSpandrels_filtered=[];
YSpandrels_filtered=[];
ZSpandrels_filtered=[];
j_spandrel=1;
% Calculate surface areas for each surface
for i = 1:num_surfaces_spandrels

    spandrel_corner1=[XSpandrels(1,i), YSpandrels(1,i), ZSpandrels(1,i)];
    spandrel_corner2=[XSpandrels(2,i), YSpandrels(2,i), ZSpandrels(2,i)];
    spandrel_corner3=[XSpandrels(3,i), YSpandrels(3,i), ZSpandrels(3,i)];
    spandrel_corner4=[XSpandrels(4,i), YSpandrels(4,i), ZSpandrels(4,i)];

    diag(i) = sqrt(sum((permute(spandrel_corner1, [1 3 2]) - permute(spandrel_corner3, [3 1 2])) .^ 2, 3));
    spandrel_vertices1(i)=sqrt(sum((permute(spandrel_corner1, [1 3 2]) - permute(spandrel_corner2, [3 1 2])) .^ 2, 3));
    spandrel_vertices2(i)=sqrt(sum((permute(spandrel_corner1, [1 3 2]) - permute(spandrel_corner4, [3 1 2])) .^ 2, 3));

    if diag(i) == spandrel_vertices1(i) || diag(i) == spandrel_vertices2(i)

    else
        j_spandrel=j_spandrel+1;
        XSpandrels_filtered(:,j_spandrel)=XSpandrels(:,i);
        YSpandrels_filtered(:,j_spandrel)=YSpandrels(:,i);
        ZSpandrels_filtered(:,j_spandrel)=ZSpandrels(:,i);

    end
    
end

% Preallocate surface areas array
num_surfaces_piers = size(XPiers, 2);
surfaceAreas_piers = zeros(1, num_surfaces_piers);

diag=zeros(size(surfaceAreas_piers));

XPiers_filtered=[];
YPiers_filtered=[];
ZPiers_filtered=[];
j_pier=1;

% Calculate surface areas for each surface
for i = 1:num_surfaces_piers

    pier_corner1=[XPiers(1,i), YPiers(1,i), ZPiers(1,i)];
    pier_corner2=[XPiers(2,i), YPiers(2,i), ZPiers(2,i)];
    pier_corner3=[XPiers(3,i), YPiers(3,i), ZPiers(3,i)];
    pier_corner4=[XPiers(4,i), YPiers(4,i), ZPiers(4,i)];

    diag(i) = sqrt(sum((permute(pier_corner1, [1 3 2]) - permute(pier_corner3, [3 1 2])) .^ 2, 3));
    pier_vertices1(i)=sqrt(sum((permute(pier_corner1, [1 3 2]) - permute(pier_corner2, [3 1 2])) .^ 2, 3));
    pier_vertices2(i)=sqrt(sum((permute(pier_corner1, [1 3 2]) - permute(pier_corner4, [3 1 2])) .^ 2, 3));

    if diag(i) == pier_vertices1(i) || diag(i) == pier_vertices2(i)

    else
        j_pier=j_pier+1;
        XPiers_filtered(:,j_pier)=XPiers(:,i);
        YPiers_filtered(:,j_pier)=YPiers(:,i);
        ZPiers_filtered(:,j_pier)=ZPiers(:,i);

    end
    
end
