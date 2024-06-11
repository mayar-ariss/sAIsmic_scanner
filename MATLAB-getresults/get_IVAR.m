function [ WCR,DAM ] = get_IVAR( element,gp,file,name )

% Gets the stresses and strains of input "element" at the selected 'gp' (gauss points) 
% for each step, from the output file name#file.num

% - The output SIG is a cell and contains one matrix at SIG{gpi} for each gp.
%   SIG{gpi}(steps,3) contains one row for each step. 
%   Each row contains the stresses [sign,tau_x,tau_y] of gpi at the 
%   corresponding step of the analysis
% - The output EPS is a cell and contains one matrix at EPS{gpi} for each gp.
%   EPS{gpi}(steps,3) contains one row for each step. 
%   Each row contains the stresses [eps_n,eps_x,eps_y] of gpi at the 
%   corresponding step of the analysis

if nargin < 4
    name='struc';
end

% Put elements and gp in order (aukson arithmos)
element=sort(element);
gp=sort(gp);

% total number of steps
nst=533;  %steps(file,name);

% Initiate SIG, EPS
WCR=cell(length(gp),length(element));
DAM=cell(length(gp),length(element));

% Open num file
if file>0
    fid=fopen(sprintf('../%01s#%03d.num',name,file),'r');
else
    fid=fopen(sprintf('../%01s.num',name),'r');
end

% Read line by line until oneline=-1 (end of file)
oneline=fgets(fid);
step=0;
key1='#ie82s1';
key2='#ie82s2';

while ischar(oneline) && step<=nst
    if sum(regexp(oneline,key1))>0
        step=step+1;
        if step>nst
            break
        end
        oneline=fgets(fid);
        oneline=fgets(fid);
        % Read every line of eps output until you find all the elements
        % in 'element'
        ie=1;
        nvar=fgets(fid);
        nvar=str2num(nvar);
        while ie<=length(element)
            if nvar(1)==element(ie);
                ig=1;
                while ig<=length(gp)
                    if nvar(2)==gp(ig)
                        DAM{ig,ie}(step,:)=nvar([9:11]);
                        ig=ig+1;
                    end
                    oneline=fgets(fid);
                    nvar=str2num(oneline);
                end
                ie=ie+1;
            else
                nvar=fgets(fid);
                nvar=str2num(nvar);
            end
        end
    end
    if sum(regexp(oneline,key2))>0
        oneline=fgets(fid);
        oneline=fgets(fid);
        % Read every line of eps output until you find all the elements
        % in 'element'
        ie=1;
        nvar=fgets(fid);
        nvar=str2num(nvar);
        while ie<=length(element)
            if nvar(1)==element(ie);
                ig=1;
                while ig<=length(gp)
                    if nvar(2)==gp(ig)
                        WCR{ig,ie}(step,:)=nvar([9:11]);
                        ig=ig+1;
                    end
                    oneline=fgets(fid);
                    nvar=str2num(oneline);
                end
                ie=ie+1;
            else
                nvar=fgets(fid);
                nvar=str2num(nvar);
            end
        end
    end
    oneline=fgets(fid);
end
fclose(fid);

% Add first zero step
for ig=1:length(gp)
    for ie=1:length(element)
        WCR{ig,ie}=vertcat([0,0,0],WCR{ig,ie});
        DAM{ig,ie}=vertcat([0,0,0],DAM{ig,ie});
    end
end

