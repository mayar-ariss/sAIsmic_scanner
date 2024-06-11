function [ EN,NC ] = read_elmcon( nbe,grpname,fname,file )
% Gets the element names and the nodal connectivity for the macroelements
% ! It assumes that: !
%        - the element connectivity is declared in one line
%        - the element name is a number
%        - the group name is declared between the element name and the
%          element connectivity

% Input:
% nbe = number of elements
% grpname = group names assigned to the macroelements we want to obtain
% fname = name of the input file name.dat or name#00i.dat
% file = number of partition (0 for parent file)


% If file number is not declares, we assume that it is the parent
if nargin < 4
    file=0;
end

% Define size of output matrices
EN=zeros(nbe,1);
NC=zeros(nbe,8);

% Open file
if file>0
    fid=fopen(sprintf('../%01s#%03d.dat',fname,file),'r');
else
    fid=fopen(sprintf('../%01s.dat',fname),'r');
end

% Initiate variables
iel=0;
fin=false;

% Start reading
oneline=fgets(fid);

while ischar(oneline)
    oneline=fgets(fid);
    
    % Find header structural.nodal
    if ~isempty(strfind(oneline,'connectivity'))
        while ~fin
            % Read each line under header. If it contains the
            % group name we are looking for store element
            oneline=fgets(fid);
            is_gp=0;
            for gp=1:length(grpname)
                if ~isempty(strfind(oneline,grpname{gp})) && ...
                        isempty(strfind(oneline,'#'))
                    is_gp=gp;
                    pos=strfind(oneline,grpname{gp});
                    break
                end
            end
            if is_gp>0
                iel=iel+1;
                oneline(pos:pos+length(grpname{is_gp}))=[];
                oneline=str2num(oneline);
                EN(iel)=oneline(1);
                NC(iel,:)=oneline(2:end);

                % When we have stored all the elements finish loop
                if iel==nbe
                    fin=true;
                end
            end
        end
        % Stop reading file
        break
    end
end
% Close the file
fclose(fid);
end
