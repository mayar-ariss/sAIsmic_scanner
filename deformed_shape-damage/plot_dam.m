% Example script for plotting deformed shape of URM macroelements in 3D
% Created 21/Nov/2016, E. Minga
% Edited Oct/2023, L. Bomben
% Edited Mar/2024, M. Ariss

%%
close all
clc 
clear all
addpath('../MATLAB-plot_adaptic'); %careful, / for macOS and \ for Windows
addpath('../MATLAB-getresults'); %careful, / for macOS and \ for Windows

%% INPUT
%%---------------------------------------------------------------------------
% fname = file name (fname.dat, fname.num or fname#00i.dat, fname#00i.num)
fname = 'W1X';
% file(:) = partition number where macroelements are declared (can be an array)
file = [0];     % [file 1,file 2]
% nbn(:) = number of nodes declared in each read file (can be an array)
nbn = [90];  % [file 1,file 2] W1X 90
% nbe(:) = total number of macroelements that we get from each file (can be an array)
nbe = [16];  % [file 1,file 2] W1X 16
% nbgp = number of groups assigned to macroelements
nbgp = 2;
% grpname(:) = group name assigned to macroelements that we want to plot (can be an array)
grpname = cell(nbgp,1);
grpname{1} = 'macro1'; 
grpname{2} = 'macro2';

% Thickness of plotted walls:
thick=250;
%vertical axis
Vert=1; %put 1 if global vertical axis==z, otherwise any other number

%damage to be plotted
%0: plotting; 1: not plotting
plot_dams=0; %in-plane
plot_damo=1; %out-of-plane

% Plot options
% plot damage in undeformed shape
plot_undef=true;
% Scale of deformed shape
plotscale=10;
% element names for damage plotting
DEN = cell(length(file),1);
for i=1
%     DEN{i}=[10101,10201,10301,10401,...
%         10108,10208,10308,10408];
    DEN{i} = 'all';
end
%%-------------------------------------------------------------------------

%% Read dat file
% NN: nodal names
% XYZ: nodal coordinated - corresponding to NN line by line
% EN: element names
% CON: element nodal connectivity - corresponding to EN line by line
% ACON: element adof connectivity - corresponding to EN line by line
[NN,XYZ,EN,CON, ACON, SPRP,DIM_EL]=read_dat_file(fname,file,nbn,nbe,grpname);

nst=steps(0,fname); %nst: final step
n=1; % n: number of plots
i=nst; % i: steps to be plotted (vector if n>1)

for j=1:n
	[DAMS,DAMO] = get_sprdam(SPRP,EN,DEN,file,fname,DIM_EL,i(j),plot_dams,plot_damo);
	%% Plot damage in undeformed shape
	if plot_undef
		% to type DAMS if in-plane damage, DAMO if out-of-plane damage
		dam_undeformed(NN,XYZ,CON,EN,DEN,DAMS,thick,Vert);
		% x,y,z plot limits
		xlim([-200 1500]);
		ylim([0 2100]);
		% printing
		%name=sprintf('Spring_damage_%1.0f_step_%1.0f.jpg',j,i(j))
		%print(gcf,fullfile(name),'-dpng','-r300')
		%close
	else
	%% Plot damage in deformed shape
	% Read num file
		DispRot=get_DispRot(file,fname);
		CADOF=get_CAdofs(file,fname,nbe,plotscale);
	
	% Create matrix with current positions and rotations
		CXYZR=get_Cxyzr(NN,XYZ,DispRot,plotscale);

	% Plot damage in deformed shape
		beamrot=cell(length(file),1);
		beamrot{1}=true; beamrot{2}=true;
		% shear springs
		dam_deformed(NN,XYZ,CXYZR,CON,beamrot,CADOF,EN,DEN,DAMS,thick,length(CADOF{1}));
		% out-of-plane springs
		%dam_deformed(NN,XYZ,CXYZR,CON,beamrot,CADOF,EN,DEN,DAMO,thick,length(CADOF{1}));
	end
end
