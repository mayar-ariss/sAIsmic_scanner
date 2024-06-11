function [ adof ] = read_cadofs( file,name,nbe,nst )

% Obtain 8 additional dofs of each macroelement of the partition at each
% step
% Input:
% - file = number of partition (0 for parent)
% - name = name of dat file
% - nbe = number of macroelements in the partition
% Output:
% - adof(nelm,8,nstep) = matrix containing the 8 adof values of each element  
%                        at each step
%                        nbe = number of macroelements in the partition
%                        8 = number of adofs per element
%                        nstep = number of steps performed in the analysis

%nst=steps(file,name);

adof=zeros(nbe,8,nst);

if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

oneline=fgets(fid);
istep=0;
while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the displacement of the nodes in 'node' (disp{node}) at the current step (ilf)
    if ~isempty(strfind(oneline,'#ie25e'))
        istep=istep+1
        oneline=fgets(fid); oneline=fgets(fid); oneline=fgets(fid);
        ielm=0;
        % Read every line of this elm output section 
        while ielm<nbe
            ielm=ielm+1;
            numline=str2num(oneline);
            if length(numline)==9
                adof(ielm,1:4,istep)=numline(end-3:end);
            else
                adof(ielm,:,istep)=numline(6:end);
            end
            oneline=fgets(fid);
        end
        if istep==nst
            fclose(fid);
            break
        end
    end
end

end

