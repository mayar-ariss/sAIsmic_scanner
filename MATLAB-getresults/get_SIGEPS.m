function [EPS,SIG] = get_SIGEPS( element,gp,file,name,nst )

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

if nst<500
    % Initiate SIG, EPS
    SIG=cell(length(gp),length(element));
    EPS=cell(length(gp),length(element));

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
            step=step+1
            if step>nst
                break
            end
            oneline=fgets(fid);
            oneline=fgets(fid);
            % Read every line of eps output until you find all the elements
            % in 'element'
            nvar=fgets(fid);
            nvar=str2num(nvar);
            ie=1;
            while ie<=length(element)
                if nvar(1)==element(ie);
                    ig=1;
                    while ig<=length(gp)
                        if nvar(2)==gp(ig)
                            EPS{ig,ie}(step,:)=nvar([5,7,8]);
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
                            SIG{ig,ie}(step,:)=nvar([5,7,8]);
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
            SIG{ig,ie}=vertcat([0,0,0],SIG{ig,ie});
            EPS{ig,ie}=vertcat([0,0,0],EPS{ig,ie});
        end
    end

else
    SIG=[];
    EPS=[];
    % Create output files and add first zero step
    fid=fopen('eps.txt','w');
    fprintf(fid, '%6.10f %6.10f %6.10f\n', [0,0,0]);
    fclose(fid);
    fid=fopen('sig.txt','w');
    fprintf(fid, '%6.10f %6.10f %6.10f\n', [0,0,0]);
    fclose(fid);
    
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
            step=step+1
            if step>nst
                break
            end
            oneline=fgets(fid);
            oneline=fgets(fid);
            % Read every line of eps output until you find element(1)
            ie=1;
            ig=1;
            fin=false;
            while ~fin
                nvar=fgets(fid);
                nvar=str2num(nvar);
                if nvar(1)==element(ie) && nvar(2)==gp(ig)
                    fid2 = fopen('eps.txt', 'a+');
                    fprintf(fid2, '%6.10f %6.10f %6.10f\n', nvar([5,7,8]));
                    fclose(fid2);
                    oneline=fgets(fid);
                    fin=true;
                end
            end
        end
        if sum(regexp(oneline,key1))>0
            step=step+1
            if step>nst
                break
            end
            oneline=fgets(fid);
            oneline=fgets(fid);
            % Read every line of eps output until you find element(1)
            ie=1;
            ig=1;
            fin=false;
            while ~fin
                nvar=fgets(fid);
                nvar=str2num(nvar);
                if nvar(1)==element(ie) && nvar(2)==gp(ig)
                    fid2 = fopen('sig.txt', 'a+');
                    fprintf(fid2, '%6.10f %6.10f %6.10f\n', nvar([5,7,8]));
                    fclose(fid2);
                    oneline=fgets(fid);
                    fin=true;
                end
            end
        end
    end
    fclose(fid);
end
