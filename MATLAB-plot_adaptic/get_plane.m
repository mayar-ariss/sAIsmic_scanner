function [ plane,T3D1,T3D2 ] = get_plane( xyz0 )
% Obtain the plane of the element from the initial coordinates of the corner nodes

x=[1,0,0];
y=[0,1,0];
z=[0,0,1];

s12=xyz0(2,:)-xyz0(1,:);
s12=round(s12./norm(s12));
s23=xyz0(3,:)-xyz0(2,:);
s23=round(s23./norm(s23));
sout=cross(s12,s23);

if s12==x
    plane(1)=1;
elseif s12==y
    plane(1)=2;
elseif s12==z
    plane(1)=3;
else
    plane(1)=1;
end
if s23==x
    plane(2)=1;
elseif s23==y
    plane(2)=2;
elseif s23==z
    plane(2)=3;
else
    plane(2)=2;
end

T3D1=[s12;s23;sout];
T3D2=[s23;-s12;sout];

end

