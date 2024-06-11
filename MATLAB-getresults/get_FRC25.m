function [FRC] = get_FRC25( element,side,file,name,nst )

% Gets the stresses and strains of of one macorelement specified in 'element'
% of one specified interface 'int' (1:4) at the selected 'gp' (gauss points, may be an array) 
% for each step, from the output file name#file.num

% Input:
% - element: an integer containing the element number
% - side: an integer from 1 to 4 specifying one side/interface
% - file: the number of the partition
% - name: the name of the dat files
% - nst: the number of steps of the analysis 
% Output:
% - FRC{nst,12} is a matrix containing the 12 forces of the given side 
%   at each step of the analysis

switch side
    case(1)
        key1='#ie25a';
    case(2)
        key1='#ie25b';
    case(3)
        key1='#ie25c';
    case(4)
        key1='#ie25d';
end
        
% Initiate SIG, EPS
FRC=zeros(nst,12);

% Open num file
if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

% Read line by line until oneline=-1 (end of file)
step=0;
fin=false;
    
    while ~fin
        oneline=fgets(fid);
        if ~isempty(strfind(oneline,key1))
            step=step+1
            oneline=fgets(fid);
            oneline=fgets(fid); oneline=fgets(fid);
            nvar=str2num(oneline);
            % Read every line of eps output until you find all the elements
            % in 'element'
            while ~fin
                if nvar(1)==element;
                    if length(nvar)==7
                        FRC(step,1:6)=nvar(2:end);
                    else
                        FRC(step,:)=nvar(2:end);
                    end
                    break
                else
                    oneline=fgets(fid);
                    nvar=str2num(oneline);
                end
            end
        end
        if step==nst
            fin=true;
        end
    end
fclose(fid);


end
