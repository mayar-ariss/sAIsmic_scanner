function [ MAS ] = matprop_mas( file,name )
% Gets the displacements of the input "node" for each step, from the output 
% file struc#file.num

% The output disp is a matrix and contains the displacements in x,y and z 
% of the node at each step

if nargin < 2
    name='struc';
end

if file>0
    fid=fopen(sprintf('../%01s#%03d.dat',name,file),'r');
else
    fid=fopen(sprintf('../%01s.dat',name),'r');
end

oneline=fgets(fid);

while ischar(oneline)
    oneline=fgets(fid);
    
    % Get the interface material properties
    if length(oneline)>4 && strcmp(oneline(1:5),' ibed')
        oneline=oneline(26:end-2);
        ibed=str2num(oneline)';
        for line=1:7
            oneline=fgets(fid);
            oneline=oneline(26:end-2);
            ibed=[ibed;str2num(oneline)'];
        end
    end
    if length(oneline)>4 && strcmp(oneline(1:6),' ihead')
        oneline=oneline(26:end-2);
        ihead=str2num(oneline)';
        for line=1:7
            oneline=fgets(fid);
            oneline=oneline(26:end-2);
            ihead=[ihead;str2num(oneline)'];
        end
    end
    if length(oneline)>4 && strcmp(oneline(1:6),' iface')
        oneline=oneline(26:end-2);
        iface=str2num(oneline)';
        for line=1:7
            oneline=fgets(fid);
            oneline=oneline(26:end-2);
            iface=[iface;str2num(oneline)'];
        end
    end
    if length(oneline)>4 && strcmp(oneline(1:7),' sbrick')
        oneline=oneline(26:end);
        sbrick=str2num(oneline)';
        break
    end
    if length(oneline)>4 && strcmp(oneline(1:4),'mat1')
        oneline=oneline(12:44);
        mat1=str2num(oneline)';
        for line=1:3
            oneline=fgets(fid);
            oneline=oneline(12:44);
            mat1=[mat1;str2num(oneline)'];
        end
        break
    end
end
fclose(fid);

if exist('ibed')
    MAS=cell(size(ibed,1)+1,5);
    MAS{2,1}='Knormal';
    MAS{3,1}='Ktang1';
    MAS{4,1}='Ktang2';
    MAS{5,1}='COH_initial';
    MAS{6,1}='COH_res';
    MAS{7,1}='TENS_strength';
    MAS{8,1}='TENS_str_res';
    MAS{9,1}='tanf_initial';
    MAS{10,1}='tanf_res';
    MAS{11,1}='parab/twocurve';
    MAS{12,1}='damage';
    MAS{13,1}='COH_initial';
    MAS{14,1}='COH_res';
    MAS{15,1}='TENS_strength';
    MAS{16,1}='TENS_str_res';
    MAS{17,1}='tanf_initial';
    MAS{18,1}='tanf_res';
    MAS{19,1}='parab/twocurve';
    MAS{19,1}='COMPR_strength';
    MAS{20,1}='COMPR_str_res';
    MAS{end-7,1}='GF1';
    MAS{end-6,1}='GF2';
    MAS{end-5,1}='GF3';

    MAS{1,2}='ibed';
    MAS{1,3}='ihead';
    MAS{1,4}='iface';
    MAS{1,5}='sbrick';
    for i=1:size(ibed,1)
        MAS{i+1,2}=ibed(i);
    end
    for i=1:size(ihead,1)
        MAS{i+1,3}=ihead(i);
    end
    for i=1:size(iface,1)
        MAS{i+1,4}=iface(i);
    end
    for i=1:size(sbrick,1)
        MAS{i+1,5}=sbrick(i);
    end
elseif exist('mat1')
    MAS=cell(size(mat1,1)+1,2);
    MAS{2,1}='Knormal';
    MAS{3,1}='Ktang1';
    MAS{4,1}='Ktang2';
    MAS{5,1}='FT';
    MAS{6,1}='TANF';
    MAS{7,1}='CO';
    MAS{8,1}='TANFG';
    MAS{9,1}='FC';
    MAS{10,1}='Gf1';
    MAS{11,1}='Gf2';
    MAS{12,1}='Gf3';
    MAS{13,1}='TOLL';
    MAS{14,1}='Type_D';
    MAS{15,1}='W2_DT';
    MAS{1,2}='mat1';
    for i=1:size(mat1,1)
        MAS{i+1,2}=mat1(i);
    end
end
