function [ R ] = get_react( node,file,name,nst,type )
% Gets the total reaction force R at the nodes in the list % "node" for 
% each step, from the output file struc#file.num

% The output R is a matrix containing the total reaction force in x,y and z from all the
% nodes in "node" at each step

if nargin < 5
    type='cont';
end

% total number of steps
if strcmp(type,'cont')
    ndof=3;
elseif strcmp(type,'full')
    ndof=6;
end
R=zeros(nst,ndof);

% open file
if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

oneline=fgets(fid);
step=1;
while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the reaction forces of the nodes in 'node' (disp{node}) 
    % at the current step (ilf)
    if length(oneline)>3 && strcmp(oneline(1:4),'#in1')
        inod=1;
        oneline=fgets(fid); oneline=fgets(fid);
        fin=false;
        % Read every line of reaction forces until you find all the
        % nodes in 'node'
        while ~fin
            nreact=fgets(fid);
            nreact(strfind(nreact, '- ')) = '0';
            if strfind(nreact, 'Infinity')>0
                pos=strfind(nreact, 'Infinity') ;
                for inf=1:length(pos)
                    nreact(pos(inf):pos(inf)+7) = '0       ';
                end
            end
            nnreact=str2num(nreact);
            if any(nnreact(1)==node)
                R(step,:)=R(step,:)+nnreact(2:1+ndof);
                inod=inod+1;
                if inod>length(node)
                    fin=true;
                end
            end
        end
        step=step+1;
        if step>nst
            break
        end
    end
    
end
fclose(fid);

end

