function [ disp ] = get_multdisp( node,file )
% Gets the displacements of the nodes in the list
% "node" for each step, from the output file struc#file.num

% The output disp is a cell array. Each cell corresponds to one of the
% nodes in "node" and contains the displacements in x,y and z of the node
% at each step

if file>0
    fid=fopen(sprintf('../struc#%03d.num',file),'r');
else
    fid=fopen('../struc.num','r');
end

oneline=fgets(fid);
step=0;
while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the displacement of the nodes in 'node' (disp{node}) at the current step (ilf)
    if length(oneline)>3 && strcmp(oneline(1:4),'#in2')
        step=step+1;
        inod=1;
        oneline=fgets(fid); oneline=fgets(fid);
        fin=false;
        % Read every line of nodal displacement until you find all the
        % nodes in 'node'
        while ~fin
            ndisp=fgets(fid);
            ndisp=str2num(ndisp);
            if ndisp(1)==node(inod);
                disp{inod}(step,:)=ndisp(2:end);
                inod=inod+1;
                if inod>size(node,1)
                    fin=true;
                end
            end
        end
    end
    
end
fclose(fid);

end
