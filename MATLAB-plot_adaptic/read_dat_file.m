function [ NN,XYZ,EN,CON, ACON, SPRP ,DIM_EL] = read_dat_file( fname,file,nbn,nbe,grpname)
%Obtains the structural nodal data and the element connectivity data from
%the adaptic input file

%   INPUT:
% fname = file name (fname.dat, fname.num or fname#00i.dat, fname#00i.num)
% file(:) = partition number where macroelements are declared (can be an array)
% nbn(:) = number of nodes declared in the file (can be an array)
% nbe(:) = total number of macroelements that we get from each file (can be an array)
% grpname(:) = group name assigned to macroelements that we want to plot 
%           (can be an array)

NN=cell(length(file),1);
XYZ=cell(length(file),1);
EN=cell(length(file),1);
CON=cell(length(file),1);
ACON=cell(length(file),1);
SPRP=cell(length(file),1);

for f=1:length(file)
    SPRP{f}=read_prop(fname,file(f));
    [NN{f},XYZ{f}]=read_coord(nbn(f),fname,file(f));
    [EN{f},CON{f}]=read_elmcon(nbe(f),grpname,fname,file(f));
	
end
for f=1:length(file)
    ACON{f}=zeros(length(EN{f}),8);
    for ll=1:length(EN{f})
        ACON{f}(ll,:)=8*ll-7:8*ll;
    end
end

for f=1:length(file)
CONN_MOD{f}=CON{f}(:,[1 3 5]);
DIMENSIONS{f}=[NN{f},XYZ{f}];
[m,n]=size(CONN_MOD{f});
DIM_EL{f}=zeros(m,2);

for i=1:m
	x1=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,1)),2);
	x2=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,2)),2);
	y1=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,3)),3);
	y2=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,1)),3);
	DIM_EL{f}(i,1)=abs(x2-x1);
	DIM_EL{f}(i,2)=abs(y2-y1);
	if DIM_EL{f}(i,1)==0 || DIM_EL{f}(i,2)==0 %direction x-z
		x1=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,1)),2);
		x2=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,2)),2);
		y1=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,3)),4);
		y2=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,1)),4);
		DIM_EL{f}(i,1)=abs(x2-x1);
		DIM_EL{f}(i,2)=abs(y2-y1);
		if DIM_EL{f}(i,1)==0 || DIM_EL{f}(i,2)==0 %direction y-z
			x1=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,1)),3);
			x2=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,2)),3);
			y1=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,3)),4);
			y2=DIMENSIONS{f}(find(DIMENSIONS{f}(:,1)==CONN_MOD{f}(i,1)),4);
			DIM_EL{f}(i,1)=abs(x2-x1);
			DIM_EL{f}(i,2)=abs(y2-y1);
		end
	end
end
end
end

