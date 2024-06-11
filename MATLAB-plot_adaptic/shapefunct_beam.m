function [ N ] = shapefunct_beam( xi,L )
% Gives matrix of shape functions for a beam in 3D with 5 dofs per node
% (3 translations, 2 rotations), given the natural coordinate xi 

% Input:
% xi: natural coordinate of beam - ATTENTION: xi=-1:1
% L: length of beam
% Output:
% N(3,10): matrix containing shape functions

% shape functions
% Nx1=0.5*(1-xi); 
% Nx2=0.5*(1+xi); 
% Ny1=(1-xi).^2.*(2+xi)./4; 
% Ny2=(1+xi).^2.*(2-xi)./4; 
% Nz1=Ny1; 
% Nz2=Ny2; 
% Nry1=L*(1-xi).^2.*(1+xi)./8;
% Nry2=L*(1+xi).^2.*(xi-1)./8; 
% Nrz1=Nry1; 
% Nrz2=Nry2;

N=zeros(3,10);
N(1,1)=0.5*(1-xi);
N(1,6)=0.5*(1+xi);
N(2,2)=(1-xi)^2*(2+xi)/4;
N(2,7)=(1+xi)^2*(2-xi)/4;
N(3,3)=N(2,2);
N(3,8)=N(2,7);
N(3,4)=L*(1-xi).^2.*(1+xi)./8;
N(3,9)=L*(1+xi).^2.*(xi-1)./8;
N(2,5)=N(3,4);
N(2,10)=N(3,9);

end

