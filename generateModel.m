% Mayar Ariss 25 Mar 2024
function [successOutput, XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face] = drawModelTremuri( model, varargin )

%[successOutput, XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered, X_face, Y_face]=generateModel(readTremuriInput('house1.txt'), 'ColorPiers', [0.2 0.2 1.0], 'ColorSpandrels', [0.2 1. .2], 'styleNodes', 'sk')
addpath('/Users/mayar/Desktop/ICL/4Y/Thesis/Codes/FEM_buildings/src/utils/hex_and_rgb_v1.1.1'); %careful, / for macOS and \ for Windows

p = inputParser;

addRequired(p,'model');
addOptional(p,'results',  [],     @isstruct);
addParameter(p,'style', 'wireframe',     @ischar);
addParameter(p,'fig', gcf,     @ishandle);

wallsToPlot = [];
for k=1:length(model.wall)
    if ~isempty(model.wall(k).x0)
        wallsToPlot = [wallsToPlot; k];
    end
end
addParameter(p,'walls', wallsToPlot,     @isnumeric);
addParameter(p,'colorPiers',     [1 1 1]*0.5);
addParameter(p,'colorSpandrels', [1 1 1]*0.5);
addParameter(p,'colorNodes',     [1 1 1]*0.5,     @isnumeric);
addParameter(p,'colorENodes',    [1 1 1]*0.0,     @isnumeric);
addParameter(p,'colorEdges',     [1 1 1]*1.0,     @isnumeric);
addParameter(p,'colorBeams',     [1 1 1]*0.5,     @isnumeric);
addParameter(p,'colorFloors',    [1 1 1]*0.7,     @isnumeric);
addParameter(p,'styleNodes', 'none',     @ischar);  % none, markers (with code)
addParameter(p,'sizeNodes', 1,     @isnumeric); 
addParameter(p,'styleBeams', 'wireframe',     @ischar);
addParameter(p,'deformed', 0,     @isnumeric);
addParameter(p,'step', 0,     @isnumeric);
addParameter(p,'elements', 0,     @isnumeric);

parse(p,model,varargin{:});

results     = p.Results.results;
style       = p.Results.style;
fig         = p.Results.fig;
wallsToPlot = p.Results.walls;
colorPiers  = p.Results.colorPiers;
colorSpandrels = p.Results.colorSpandrels;
colorNodes  = p.Results.colorNodes;
colorENodes = p.Results.colorENodes;
colorEdges  = p.Results.colorEdges;
colorBeams  = p.Results.colorBeams;
colorFloors = p.Results.colorFloors;
styleNodes  = p.Results.styleNodes;
sizeNodes   = p.Results.sizeNodes;
styleBeams  = p.Results.styleBeams;
deformed    = p.Results.deformed;
stepToPlot  = p.Results.step;
elementList = p.Results.elements;

if elementList == 0
    elementList = 1:length(model.element);
end



%%old code

%% draw
set(0, 'CurrentFigure', fig)
%axis equal

%initialise node vectors
X_vec = zeros(length(model.node2d), 1)*NaN;
Y_vec = zeros(length(model.node2d), 1)*NaN;
Z_vec = zeros(length(model.node2d), 1)*NaN;

% draw all nodes
XNodes = []; YNodes = []; ZNodes = [];
for k=1:length(model.node2d)
    if ~isempty(model.node2d(k).wall) && sum(model.node2d(k).wall == wallsToPlot)~=0
        % x and y coordinate
        if stepToPlot==0 % not deformed
            X_vec(k) = model.node2d(k).x;
            Y_vec(k) = model.node2d(k).y;
            Z_vec(k) = model.node2d(k).z;
          
        else  % deformed
            % check if we have the true hor. displacement of the node
            dirVec = [cos(model.wall(model.node2d(k).wall).angle);
                sin(model.wall(model.node2d(k).wall).angle)];
            % out of plane displacement (from 3d nodes)
            dist1 =       [model.node2d(k).x - model.node3d(model.node2d(k).repartition(1)).x;
                model.node2d(k).y - model.node3d(model.node2d(k).repartition(1)).y;
                model.node2d(k).z - model.node3d(model.node2d(k).repartition(1)).z ];
            dist2 =       [model.node2d(k).x - model.node3d(model.node2d(k).repartition(2)).x;
                model.node2d(k).y - model.node3d(model.node2d(k).repartition(2)).y;
                model.node2d(k).z - model.node3d(model.node2d(k).repartition(2)).z ];
            if isempty(results.node3d(model.node2d(k).repartition(1)).ux)
                results.node3d(model.node2d(k).repartition(1)).ux(stepToPlot) = 0;
            end
            
            if isempty(results.node3d(model.node2d(k).repartition(1)).uy)
                results.node3d(model.node2d(k).repartition(1)).uy(stepToPlot) = 0;
            end
            
            if isempty(results.node3d(model.node2d(k).repartition(1)).uz)
                results.node3d(model.node2d(k).repartition(1)).uz(stepToPlot) = 0;
            end
            
            if isempty(results.node3d(model.node2d(k).repartition(2)).ux)
                results.node3d(model.node2d(k).repartition(2)).ux(stepToPlot) = 0;
            end
            
            if isempty(results.node3d(model.node2d(k).repartition(2)).uy)
                results.node3d(model.node2d(k).repartition(2)).uy(stepToPlot) = 0;
            end
            
            if isempty(results.node3d(model.node2d(k).repartition(2)).uz)
                results.node3d(model.node2d(k).repartition(2)).uz(stepToPlot) = 0;
            end
            
            
            d1 =          [results.node3d(model.node2d(k).repartition(1)).ux(stepToPlot);
                results.node3d(model.node2d(k).repartition(1)).uy(stepToPlot);
                results.node3d(model.node2d(k).repartition(1)).uz(stepToPlot)];
            d2 =          [results.node3d(model.node2d(k).repartition(2)).ux(stepToPlot);
                results.node3d(model.node2d(k).repartition(2)).uy(stepToPlot);
                results.node3d(model.node2d(k).repartition(2)).uz(stepToPlot)];
            
            if ~isempty(results.node2d(k).ux)
                
                d1(1:2) = d1(1:2) - dirVec'*d1(1:2) *dirVec;
                d2(1:2) = d2(1:2) - dirVec'*d2(1:2) *dirVec;
                d = d1 * norm(dist2) / ( norm(dist1)+norm(dist2) ) + d2 * norm(dist1) / ( norm(dist1)+norm(dist2) );
                
                X_vec(k) = model.node2d(k).x + deformed*(d(1) + dirVec(1)*results.node2d(k).ux(stepToPlot));
                Y_vec(k) = model.node2d(k).y + deformed*(d(2) + dirVec(2)*results.node2d(k).ux(stepToPlot));
                Z_vec(k) = model.node2d(k).z + deformed*d(3);
                
            else  % calculate it from the 3d nodes. if they are not recorded, 0
                d = d1 * norm(dist2) / ( norm(dist1)+norm(dist2) ) + d2 * norm(dist1) / ( norm(dist1)+norm(dist2) );
                X_vec(k) = model.node2d(k).x + deformed*d(1);
                Y_vec(k) = model.node2d(k).y + deformed*d(2);
                Z_vec(k) = model.node2d(k).z + deformed*d(3);
            end  
        end
        
        if model.node2d(k).type == 'R'
            nodeI = [nanmean(model.node2d(k).offsetX);
                     nanmean(model.node2d(k).offsetY);
                     Z_vec(k) + model.node2d(k).offsetZ(2)];
            nodeJ = [nanmean(model.node2d(k).offsetX);
                     nanmean(model.node2d(k).offsetY);
                     Z_vec(k) + model.node2d(k).offsetZ(1)];
            
            nWall = model.node2d(k).wall;
            theta = model.wall(nWall).angle;
            v_or = [-sin(theta); cos(theta); 0];
            
            L = abs(model.node2d(k).offsetXloc(2) -  model.node2d(k).offsetXloc(1));
            t = model.node2d(k).thickness;
            
            [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, L, t, [0;0;0], [0;0;0], style);
            XNodes = [XNodes, xx];
            YNodes = [YNodes, yy];
            ZNodes = [ZNodes, zz];
        elseif model.node2d(k).type == 'P'           
            for kPolygon=1:length(model.node2d(k).polygon)
                if ~isempty(model.node2d(k).polygon(kPolygon).rho)
                    nodeI = [nanmean(model.node2d(k).polygon(kPolygon).offsetX);
                            nanmean(model.node2d(k).polygon(kPolygon).offsetY);
                            Z_vec(k) + model.node2d(k).polygon(kPolygon).offsetZ(2)];
                    nodeJ = [nanmean(model.node2d(k).polygon(kPolygon).offsetX);
                            nanmean(model.node2d(k).polygon(kPolygon).offsetY);
                            Z_vec(k) + model.node2d(k).polygon(kPolygon).offsetZ(1)];
                    
                    nWall = model.node2d(k).wall;
                    theta = model.wall(nWall).angle;
                    v_or = [-sin(theta); cos(theta); 0];
                    
                    L = abs(model.node2d(k).polygon(kPolygon).offsetXloc(2) -  model.node2d(k).polygon(kPolygon).offsetXloc(1));
                    t = model.node2d(k).polygon(kPolygon).thickness;
                    
                    [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, L, t, [0;0;0], [0;0;0], style);
                    XNodes = [XNodes, xx];
                    YNodes = [YNodes, yy];
                    ZNodes = [ZNodes, zz];             
                end
                
            end
            
        end
    end
end

%eddy changed
Xnodes=XNodes;
Ynodes=YNodes;
Znodes=ZNodes;

if ~strcmp(styleNodes, 'none')
    plot3(X_vec, Y_vec, Z_vec, styleNodes, 'MarkerFaceColor', colorENodes, 'MarkerEdgeColor', colorENodes, 'MarkerSize', sizeNodes);
    hold on
    
    if ischar(colorNodes)
        h = patch(XNodes, YNodes, ZNodes, 'w', 'edgecolor', colorEdges);
        set(h, 'FaceColor', colorNodes);
    else
        patch(XNodes, YNodes, ZNodes, colorNodes, 'edgecolor', colorEdges);
    end
    
    
end

XNodes = []; YNodes = []; ZNodes = [];
for k=1:length(model.node3d)
    if not(isempty(model.node3d(k).wall)) && (sum(model.node3d(k).wall(1) == wallsToPlot)~=0  || sum(model.node3d(k).wall(2) == wallsToPlot)~=0)
        if stepToPlot==0 % not deformed
            XX_vec(k) = model.node3d(k).x;
            YY_vec(k) = model.node3d(k).y;
            ZZ_vec(k) = model.node3d(k).z;
        else  % deformed
            
            % check if we have the true hor. displacement of the node
            if ~isempty(results.node3d(k).ux)
                XX_vec(k) = model.node3d(k).x + deformed* results.node3d(k).ux(stepToPlot);
            else  % calculate it from the 3d nodes. if they are not recorded, undeformed
                XX_vec(k) = model.node3d(k).x;
            end
            
            if ~isempty(results.node3d(k).uy)
                YY_vec(k) = model.node3d(k).y + deformed* results.node3d(k).uy(stepToPlot);
            else  % calculate it from the 3d nodes. if they are not recorded, undeformed
                YY_vec(k) = model.node3d(k).y;
            end
            
            if ~isempty(results.node3d(k).uz)
                ZZ_vec(k) = model.node3d(k).z + deformed* results.node3d(k).uz(stepToPlot);
            else  % calculate it from the 3d nodes. if they are not recorded, undeformed
                ZZ_vec(k) = model.node3d(k).z;
            end
        end
        
        if model.node3d(k).type(1) == 'R'
            nodeI = [nanmean(model.node3d(k).offsetX1);
                     nanmean(model.node3d(k).offsetY1);
                     ZZ_vec(k) + model.node3d(k).offsetZ1(2)];
            nodeJ = [nanmean(model.node3d(k).offsetX1);
                     nanmean(model.node3d(k).offsetY1);
                     ZZ_vec(k) + model.node3d(k).offsetZ1(1)];
            
            nWall = model.node3d(k).wall(1);
            theta = model.wall(nWall).angle;
            v_or = [-sin(theta); cos(theta); 0];
            
            L = abs(model.node3d(k).offsetXloc1(2) -  model.node3d(k).offsetXloc1(1));
            t = model.node3d(k).thickness1;
            
            [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, L, t, [0;0;0], [0;0;0], style);
            XNodes = [XNodes, xx];
            YNodes = [YNodes, yy];
            ZNodes = [ZNodes, zz];
            
        elseif (model.node3d(k).type(1) == 'P')          
            for kPolygon=1:length(model.node3d(k).polygon1)
                if ~isempty(model.node3d(k).polygon1(kPolygon).rho)
                    nodeI = [nanmean(model.node3d(k).polygon1(kPolygon).offsetX);
                            nanmean(model.node3d(k).polygon1(kPolygon).offsetY);
                            ZZ_vec(k) + model.node3d(k).polygon1(kPolygon).offsetZ(2)];
                    nodeJ = [nanmean(model.node3d(k).polygon1(kPolygon).offsetX);
                            nanmean(model.node3d(k).polygon1(kPolygon).offsetY);
                            ZZ_vec(k) + model.node3d(k).polygon1(kPolygon).offsetZ(1)];
                    
                    nWall = model.node3d(k).wall(1);
                    theta = model.wall(nWall).angle;
                    v_or = [-sin(theta); cos(theta); 0];
                    
                    L = abs(model.node3d(k).polygon1(kPolygon).offsetXloc(2) -  model.node3d(k).polygon1(kPolygon).offsetXloc(1));
                    t = model.node3d(k).polygon1(kPolygon).thickness;
                    
                     [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, L, t, [0;0;0], [0;0;0], style);
                     XNodes = [XNodes, xx];
                     YNodes = [YNodes, yy];
                     ZNodes = [ZNodes, zz];       
                     
                  
                end
                
            end          
            
        end
        
        
        
         if model.node3d(k).type(2) == 'R'     
            nodeI = [nanmean(model.node3d(k).offsetX2);
                nanmean(model.node3d(k).offsetY2);
                ZZ_vec(k) + model.node3d(k).offsetZ2(2)];
            nodeJ = [nanmean(model.node3d(k).offsetX2);
                nanmean(model.node3d(k).offsetY2);
                ZZ_vec(k) + model.node3d(k).offsetZ2(1)];
            
            nWall = model.node3d(k).wall(2);
            theta = model.wall(nWall).angle;
            v_or = [-sin(theta); cos(theta); 0];
            
            L = abs(model.node3d(k).offsetXloc2(2) -  model.node3d(k).offsetXloc2(1));
            t = model.node3d(k).thickness2;
            
            [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, L, t, [0;0;0], [0;0;0], style);
            XNodes = [XNodes, xx];
            YNodes = [YNodes, yy];
            ZNodes = [ZNodes, zz];
            
         elseif model.node3d(k).type(2) == 'P'
             for kPolygon=1:length(model.node3d(k).polygon2)                
                 if ~isempty(model.node3d(k).polygon2(kPolygon).rho)
                     nodeI = [nanmean(model.node3d(k).polygon2(kPolygon).offsetX);
                              nanmean(model.node3d(k).polygon2(kPolygon).offsetY);
                              ZZ_vec(k) + model.node3d(k).polygon2(kPolygon).offsetZ(2)];
                     nodeJ = [nanmean(model.node3d(k).polygon2(kPolygon).offsetX);
                         nanmean(model.node3d(k).polygon2(kPolygon).offsetY);
                         ZZ_vec(k) + model.node3d(k).polygon2(kPolygon).offsetZ(1)];
                     
                     nWall = model.node3d(k).wall(2);
                     theta = model.wall(nWall).angle;
                     v_or = [-sin(theta); cos(theta); 0];
                     
                     L = abs(model.node3d(k).polygon2(kPolygon).offsetXloc(2) -  model.node3d(k).polygon2(kPolygon).offsetXloc(1));
                     t = model.node3d(k).polygon2(kPolygon).thickness;
                     
                     [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, L, t, [0;0;0], [0;0;0], style);
                     XNodes = [XNodes, xx];
                     YNodes = [YNodes, yy];
                     ZNodes = [ZNodes, zz];
                 end
                 
             end
            
        end
        
       
        
    else
        XX_vec(k) = NaN;
        YY_vec(k) = NaN;
        ZZ_vec(k) = NaN;
    end
end


%eddy changed
Xnodes=[XNodes Xnodes];
Ynodes=[YNodes Ynodes];
Znodes=[ZNodes Znodes];


if ~strcmp(styleNodes, 'none')
    plot3(XX_vec, YY_vec, ZZ_vec, styleNodes, 'MarkerFaceColor', colorNodes, 'MarkerEdgeColor', colorENodes, 'MarkerSize', sizeNodes);
    hold on
        if ischar(colorNodes)
        h = patch(XNodes, YNodes, ZNodes, 'w', 'edgecolor', colorEdges);
        set(h, 'FaceColor', colorNodes);
    else
        patch(XNodes, YNodes, ZNodes, colorNodes, 'edgecolor', colorEdges);
    end
end

% initialise variables
XPiers = []; YPiers = []; ZPiers = [];
XSpandrels = []; YSpandrels= []; ZSpandrels = [];

for k_el=1:length(elementList)
    k = elementList(k_el);
    if not(isempty(model.element(k).angle)) && sum(model.element(k).wall == wallsToPlot)~=0
        nWall = model.element(k).wall;
        theta = model.wall(nWall).angle;
        v_or = [-sin(theta); cos(theta); 0];
        
        if stepToPlot==0 % not deformed
            nodeI = model.element(k).nI;
            nodeJ = model.element(k).nJ;
        else
            if model.element(k).nodeI > length(model.node2d)
                model.node2d(model.element(k).nodeI).wall = [];
            end
            
            %if model.element(k).nodeI <= length(model.node3d)
             if isempty(model.node2d(model.element(k).nodeI).wall)
                nI  = [model.node3d(model.element(k).nodeI).x;
                    model.node3d(model.element(k).nodeI).y;
                    model.node3d(model.element(k).nodeI).z];
                nnI = [XX_vec(model.element(k).nodeI);
                    YY_vec(model.element(k).nodeI);
                    ZZ_vec(model.element(k).nodeI)];
                hi = model.element(k).nI - nI;
                v_or_i = v_or;
                
            else
                nI = [model.node2d(model.element(k).nodeI).x;
                    model.node2d(model.element(k).nodeI).y;
                    model.node2d(model.element(k).nodeI).z];
                nnI = [X_vec(model.element(k).nodeI);
                    Y_vec(model.element(k).nodeI);
                    Z_vec(model.element(k).nodeI)];
                hi = model.element(k).nI - nI;
                
                n1 = model.node2d(model.element(k).nodeI).repartition(1);
                n2 = model.node2d(model.element(k).nodeI).repartition(2);
                dx = XX_vec(n1) - XX_vec(n2);
                dy = YY_vec(n1) - YY_vec(n2);
                l = sqrt(dx^2+dy^2);
                dx = dx/l;
                dy = dy/l;
                v_or_i = [-dy; dx; 0];
             end
            
            if model.element(k).nodeJ > length(model.node2d)
                model.node2d(model.element(k).nodeJ).wall = [];
            end
            
            %if model.element(k).nodeI <= length(model.node3d)
             if isempty(model.node2d(model.element(k).nodeJ).wall)
                nJ = [model.node3d(model.element(k).nodeJ).x;
                    model.node3d(model.element(k).nodeJ).y;
                    model.node3d(model.element(k).nodeJ).z];
                nnJ = [XX_vec(model.element(k).nodeJ);
                    YY_vec(model.element(k).nodeJ);
                    ZZ_vec(model.element(k).nodeJ)];
                hj = model.element(k).nJ - nJ;
                v_or_j = v_or;
                
            else
                nJ = [model.node2d(model.element(k).nodeJ).x;
                    model.node2d(model.element(k).nodeJ).y;
                    model.node2d(model.element(k).nodeJ).z];
                nnJ = [X_vec(model.element(k).nodeJ);
                    Y_vec(model.element(k).nodeJ);
                    Z_vec(model.element(k).nodeJ)];
                hj = model.element(k).nJ - nJ;
                
                n1 = model.node2d(model.element(k).nodeJ).repartition(1);
                n2 = model.node2d(model.element(k).nodeJ).repartition(2);
                dx = XX_vec(n1) - XX_vec(n2);
                dy = YY_vec(n1) - YY_vec(n2);
                l = sqrt(dx^2+dy^2);
                dx = dx/l;
                dy = dy/l;
                v_or_j = [-dy; dx; 0];
            end
            
            if model.element(k).type == 0
                if v_or_i'*v_or_j>0
                    v_or = (v_or_i + v_or_j)/2;
                else
                    v_or = (-v_or_i + v_or_j)/2;
                end
            end
            
            h  = (nJ-nI);
            h = h/norm(h);
            
            hi_0 =  h * h'*hi;
            hi_90 = hi - hi_0;
            hj_0 =  h * h'*hj;
            hj_90 = hj - hj_0;
            
            hh  = (nnJ-nnI); hh = hh/norm(hh);
            
            nodeI =  nnI + norm(hi_0)*hh/norm(hh) + hi_90;
            nodeJ =  nnJ - norm(hj_0)*hh/norm(hh) + hj_90;
        end
        
        [xx yy zz] = addWallToPatch(nodeI, nodeJ, v_or, model.element(k).L, model.element(k).t, [0;0;0], [0;0;0], style);
        
        if model.element(k).type == 0
            %if strcmp(style, 'wireframe'); colorEdges = colorPiers;  end
            %draw_wall( nodeI, nodeJ, v_or, model.element(k).L, model.element(k).t, ...
            %           [0;0;0], [0;0;0], style, colorEdges, colorPiers);
            XPiers = [XPiers, xx];
            YPiers = [YPiers, yy];
            ZPiers = [ZPiers, zz];
        else
            %if strcmp(style, 'wireframe'); colorEdges = colorSpandrels;  end
            %draw_wall( nodeI, nodeJ, v_or, model.element(k).L, model.element(k).t, ...
            %           [0;0;0], [0;0;0], style, colorEdges, colorSpandrels);
            XSpandrels = [XSpandrels, xx];
            YSpandrels = [YSpandrels, yy];
            ZSpandrels = [ZSpandrels, zz];
        end
    end
end

if ischar(colorPiers)
    h = patch(XPiers, YPiers, ZPiers, 'w', 'edgecolor', colorEdges);
    set(h, 'FaceColor', colorPiers);
else
   patch(XPiers, YPiers, ZPiers, colorPiers, 'edgecolor', colorEdges);
end

if ischar(colorPiers)
    h = patch(XSpandrels, YSpandrels, ZSpandrels, 'w', 'edgecolor', colorEdges);
    set(h, 'FaceColor', colorSpandrels);
else
    patch(XSpandrels, YSpandrels, ZSpandrels, colorSpandrels, 'edgecolor', colorEdges);
end


%axis equal
%axis off


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

num_surfaces_nodes = size(Xnodes, 2);
surfaceAreas_nodes = zeros(1, num_surfaces_nodes);

node_diag=zeros(size(surfaceAreas_nodes));

XNodes_filtered=[];
YNodes_filtered=[];
ZNodes_filtered=[];
j_node=0;

% Calculate surface areas for each surface
for i = 1:num_surfaces_nodes

    node_corner1=[Xnodes(1,i), Ynodes(1,i), Znodes(1,i)];
    node_corner2=[Xnodes(2,i), Ynodes(2,i), Znodes(2,i)];
    node_corner3=[Xnodes(3,i), Ynodes(3,i), Znodes(3,i)];
    node_corner4=[Xnodes(4,i), Ynodes(4,i), Znodes(4,i)];

    node_diag(i) = sqrt(sum((permute(node_corner1, [1 3 2]) - permute(node_corner3, [3 1 2])) .^ 2, 3));
    node_vertices1(i)=sqrt(sum((permute(node_corner1, [1 3 2]) - permute(node_corner2, [3 1 2])) .^ 2, 3));
    node_vertices2(i)=sqrt(sum((permute(node_corner1, [1 3 2]) - permute(node_corner4, [3 1 2])) .^ 2, 3));

    if node_diag(i) == node_vertices1(i) || node_diag(i) == node_vertices2(i)

    else
        j_node=j_node+1;
        XNodes_filtered(:,j_node)=Xnodes(:,i);
        YNodes_filtered(:,j_node)=Ynodes(:,i);
        ZNodes_filtered(:,j_node)=Znodes(:,i);

    end
    
end


% Define matrices XNodes_filtered, YNodes_filtered, XSpandrels_filtered,
% YSpandrels_filtered, XPiers_filtered, and YPiers_filtered

% Define a function to find values occurring four times consecutively
find_consecutive_four = @(matrix) unique(matrix(all(diff(matrix) == 0, 1), :));

% Initialize variables
X_face = [];
Y_face = [];
X_matrices = {round(XNodes_filtered, 4), round(XSpandrels_filtered, 4), round(XPiers_filtered, 4)};
Y_matrices = {round(YNodes_filtered, 4), round(YSpandrels_filtered, 4), round(YPiers_filtered, 4)};

% Find X_face for each matrix
for i = 1:numel(X_matrices)
    matrix = X_matrices{i};
    [rows, cols] = size(matrix);
    for col = 1:cols
        for row = 1:rows-3
            if all(matrix(row:row+3, col) == matrix(row, col))
                X_face = [X_face; matrix(row, col)];
            end
        end
    end
end

% Find Y_face for each matrix
for i = 1:numel(Y_matrices)
    matrix = Y_matrices{i};
    [rows, cols] = size(matrix);
    for col = 1:cols
        for row = 1:rows-3
            if all(matrix(row:row+3, col) == matrix(row, col))
                Y_face = [Y_face; matrix(row, col)];
            end
        end
    end
end

% Extract unique values
X_face = unique(X_face);
Y_face = unique(Y_face);



successOutput = true;


end

