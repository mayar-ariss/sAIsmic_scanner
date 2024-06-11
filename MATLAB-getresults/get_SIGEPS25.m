function [EPS,SIG] = get_SIGEPS25( element,int,gp,file,name,nst )

% Gets the stresses and strains of of one macorelement specified in 'element'
% of one specified interface 'int' (1:4) at the selected 'gp' (gauss points, may be an array) 
% for each step, from the output file name#file.num

% Input:
% - element: an integer containing the element number
% - int: an integer from 1 to 4 specifying one interface
% - gp: an array containing the gp we want to control
% - file: the number of the partition
% - name: the name of the dat files
% - nst: the number of steps of the analysis 
% Output:
% - The output SIG is a cell and contains one matrix SIG{gpi} for each gp.
%   SIG{gpi}(steps,3) contains one row for each step. 
%   Each row contains the stresses [sign,tau_x,tau_y] of gpi at the 
%   corresponding step of the analysis
% - The output EPS is a cell and contains one matrix at EPS{gpi} for each gp.
%   EPS{gpi}(steps,3) contains one row for each step. 
%   Each row contains the stresses [eps_n,eps_x,eps_y] of gpi at the 
%   corresponding step of the analysis

key1='#ie25s';

% Initiate SIG, EPS
ngp=length(gp);
SIG=cell(ngp,1);
EPS=cell(ngp,1);
for gg=1:ngp
    SIG{gg}=zeros(nst,3);
    EPS{gg}=zeros(nst,3);
end

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
                    ig=1;
                    while ig<=length(gp)
                        if nvar(2)==int+gp(ig)*10e-5
                            EPS{ig}(step,:)=nvar(3:5);
                            SIG{ig}(step,:)=nvar(end-2:end);
                            ig=ig+1;
                        end
                        oneline=fgets(fid);
                        nvar=str2num(oneline);
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

% Add first zero step
for ig=1:length(gp)
   SIG{ig}=vertcat([0,0,0],SIG{ig});
   EPS{ig}=vertcat([0,0,0],EPS{ig});
end

end
