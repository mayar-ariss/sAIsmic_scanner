function [ acc ] = get_acc( node,file,name,nst )
% Gets the displacements of the input "node" for each step, from the output 
% file struc#file.num

% The output disp is a matrix and contains the displacements in x,y and z 
% of the node at each step
if nargin < 3
    name='struc';
end

acc=zeros(nst,3*length(node));

if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

oneline=fgets(fid);
step=0;
while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the displacement of the nodes in 'node' (disp{node}) at the current step (ilf)
    if length(oneline)>3 && strcmp(oneline(1:4),'#in4')
        step=step+1
        inod=1;
        oneline=fgets(fid); oneline=fgets(fid);
        fin=false;
        % Read every line of nodal displacement until you find all the
        % nodes in 'node'
        while ~fin
            ndisp=fgets(fid);
            ndisp=str2num(ndisp);
            if any(ndisp(1)==node)
                if length(node)>1
                    pos=(ndisp(1)==node);
                    acc(step,reshape(repmat(pos,3,1),1,[]))=ndisp(2:4);
                else
                    acc(step,:)=ndisp(2:4);
                end
                inod=inod+1;
                if inod>length(node)
                    fin=true;
                end
            end
        end
        if step==nst
            fclose(fid);
            break
        end
    end
    
end

end
