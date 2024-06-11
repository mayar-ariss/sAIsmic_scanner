% Example script for plotting deformed shape of URM macroelements in 3D
% Edited Mar/2024, M. Ariss

%%
clear; clc; close all
addpath('/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/MATLAB-plot_adaptic'); %careful, / for macOS and \ for Windows
addpath('/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/MATLAB-getresults'); %careful, / for macOS and \ for Windows
addpath('/Users/mayar/Desktop/ICL/4Y/Thesis'); %careful, / for macOS and \ for Windows

%% INPUT
%%---------------------------------------------------------------------------
% fname = file name (fname.dat, fname.num or fname#00i.dat, fname#00i.num)
fname = ['W1Y'];
% file(:) = partition number where macroelements are declared (can be an array)
file = [0];     % [file 1,file 2] short_panel:0 W1Y 
% nbn(:) = number of nodes declared in each file read (can be an array)
nbn = [178];  % [file 1,file 2] W1Y 178 lacapite 104; school 176; house1 136;
% nbe(:) = total number of macroelements that we get from each file (can be an array)
nbe = [35];  % [file 1,file 2]  W1Y 16 short_panel:15; lacapite 13; school 22; house1 17;
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
% plot undeformed shape
plot_undef=true;
% Scale of deformed shape
plotscale=1000;
%%-------------------------------------------------------------------------

%% Read dat file
% NN: nodal names
% XYZ: nodal coordinated - corresponding to NN line by line
% EN: element names
% CON: element nodal connectivity - corresponding to EN line by line
% ACON: element adof connectivity - corresponding to EN line by line
[NN,XYZ,EN,CON, ACON, SPRP,DIM_EL]=read_dat_file(fname,file,nbn,nbe,grpname);

%% Plot undeformed shape
if plot_undef
    plot_undeformed(NN,XYZ,CON,thick,Vert);
     %xlim([-1 6]);
     %ylim([-1 6]);
end

%% Read num file
[DispRot,nst]=get_DispRot(file,fname);
CADOF=get_CAdofs(file,fname,nbe,nst,plotscale);

%% Create matrix with current positions and rotations
CXYZR=get_Cxyzr(NN,XYZ,DispRot(:,:,:),plotscale);

%% Plot deformed shape
beamrot=cell(length(file),1);
beamrot{1}=false;
n=1; % number of plots
i=nst; % steps to be plotted (vector if n>1)

for j=1:n
	plot_deformed(NN,XYZ,CXYZR,CON,beamrot,CADOF,thick,i(j),Vert);
	% x,y,z plot limits
	xlim([-300 3000]);
	ylim([0 2000]);
	% printing
	%name=sprintf('Def_shape_%1.0f_passo_%1.0f.jpg',j,i(j))
	%name='Def_shape_j.jpg';
	%print(gcf,fullfile(name),'-dpng','-r300')
end