function [ dxyz ] = transform3d( NXYZ,dlocal )
% Transforms adof displacements in the local (element) coordinate system to
% displacements in the global xyz coordinate system

dxyz=zeros(4,3);

% Find plane of element
x=NXYZ(2,2:end)-NXYZ(1,2:end);
y=NXYZ(7,2:end)-NXYZ(1,2:end);
x=x/norm(x);
y=y/norm(y);

% Transform
if sum(x==[1,0,0])==3
    dxyz(:,1)=dlocal(:,1);
    if sum(y==[0,1,0])==3
        dxyz(:,2)=dlocal(:,2);
        dxyz(:,3)=dlocal(:,3);
    elseif sum(y==[0,0,1])==3
        dxyz(:,3)=dlocal(:,2);
        dxyz(:,2)=dlocal(:,3);
    end
elseif sum(x==[0,1,0])==3
    dxyz(:,2)=dlocal(:,1);
    if sum(y==[0,0,1])==3
        dxyz(:,3)=dlocal(:,2);
        dxyz(:,1)=-dlocal(:,3);
    elseif sum(y==[1,0,0])==3
        dxyz(:,3)=-dlocal(:,3);
        dxyz(:,1)=dlocal(:,2);
    end 
elseif sum(x==[0,0,1])==3
    dxyz(:,3)=dlocal(:,1);
    if sum(y==[0,1,0])==3
        dxyz(:,2)=dlocal(:,2);
        dxyz(:,1)=-dlocal(:,3);
    elseif sum(y==[1,0,0])==3
        dxyz(:,1)=dlocal(:,2);
        dxyz(:,2)=-dlocal(:,3);
    end 
end
