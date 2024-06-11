function [ mtx2 ] = matrix_part( ind,col,matrix,unique)
%Extracts the part of the matrix that has the elements of the array "ind"
%in the column col
%'unique'=true(default)/false defines how multiple identical entries in ind are treated
%         if unique = true, then the order of the entries in mtx2 is taken from
%                     their order in matrix
%         if unique = false, then the order of the entries in mtx2 is taken from
%                     their order in ind

% for example:
% unique=true
% matrix = [1,0,0,0; 2,0,0,1; 3,1,0,0; 4,0,1,0]  %(matrix of nodal names and 3d coordinates)
% ind = [3,1];  %(I am interested in node numbers 1,3)
% col = 1; %(the column containing my identifying characteristic is column 1)
% mtx2 = [1,0,0,0; 3,0,0]

% unique=false (implies order='ind')
% matrix = [1,0,0,0; 2,0,0,1; 3,1,0,0; 4,0,1,0]  %(matrix of nodal names and 3d coordinates)
% ind = [3,1,3];  %(I am interested in node numbers 1,3)
% col = 1; %(the column containing my identifying characteristic is column 1)
% mtx2 = [3,1,0,0; 1,0,0,0; 3,1,0,0]

if nargin <4
    unique=true;
end
    
if unique
    ir = ismember(matrix(:,col),ind);
    mtx2 = matrix(ir,:);
else
    [~,ir]=ismember(ind,matrix(:,col));
    
    mtx2 = matrix(ir,:);
end

end

