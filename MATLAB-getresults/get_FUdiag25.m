function [FUD] = get_FUdiag25( element,file,name,nst )

% Gets the diagonal force-displacement of one macorelement specified in 'element'
% for each step, from the output file name#file.num

% Input:
% - element: an integer containing the element number
% - file: the number of the partition
% - name: the name of the dat files
% - nst: the number of steps of the analysis 
% Output:
% - FRC{nst,12} is a matrix containing the 12 forces of the given side 
%   at each step of the analysis

key1='#ie25e';
        
% Initiate SIG, EPS
FUD=cell(length(element),1);

% Open num file
if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

% Read line by line until oneline=-1 (end of file)
step=0;
fin=false;

% commented lines were aimed to plot F-d cycles of diagonal spring 
%name = 'F-d_diagonal_spring.txt';
%fileID = fopen(sprintf('dat_files/%01s.dat',name),'w');
%x = [0:.1:1];
%A = [x; exp(x);(-x)];

%fileID = fopen(name,'w');
%fprintf(fileID,'%8s %4s %12s %12s\n','Elm.name','step','Fd1','ud1');

    while ~fin
        oneline=fgets(fid);
        if ~isempty(strfind(oneline,key1))
            step=step+1
            oneline=fgets(fid);
            oneline=fgets(fid); oneline=fgets(fid);
            ne=0;
            nvar=str2num(oneline);
            % Read every line of eps output until you find all the elements
            % in 'element'
            while ~fin
                [~,ee]=ismember(nvar(1),element);
                if ee~=0
                    ne=ne+1;
                    FUD{ee}(step,:)=nvar(2:5);
					%fprintf(fileID,'%8.0f %4.0f %12.8f %12.8f\n',nvar(1),step,nvar(2:3));
					if ne==length(element)
                        break
                    else
                        oneline=fgets(fid);
                        nvar=str2num(oneline);
                    end
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
%fclose(fileID);

end
