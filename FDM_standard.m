% 2D elasticity FDM ---
% 
% Mohsen Goudarzi (2020 - Utrecht) 
% ---------------------------------------

addpath ('/home/mohsen/Downloads/XFEMMatlabcode')
addpath ('/home/mohsen/Desktop/MHSN/MatlabTools/p1top2')
clear all ; 
close all ; 
clc
tic;     


% Dimension of the domain 
D = 4 ;
L = 4 ;

% Material properties
E  = 1e2 ;
nu = 0.2 ;
stressState='PLANE_STRAIN';


% Number of nodes along two directions
nnx = 10 ;
nny = 10 ;

% Four corner points
pt1 = [0 0] ; pt2 = [D 0] ; pt3 = [D L] ; pt4 = [0 L] ;

elemType = 'Q4' ;
[node,element] = meshRectangularRegion(...
    pt1, pt2, pt3, pt4, nnx,nny,elemType);

figure
hold on
plot(node(:,1),node(:,2),'rsq')
axis equal

dX = norm(node(2,1)-node(1,1));
dY = dX ; 

%% uy DOFs
    node_y = [] ; 
for ij = 1 : size(node,1)
    if rem(ij,nnx) == 1
        node_new = [-dX/2 node(ij,2)] ;
        node_y = [ node_y ; [-dX/2 node(ij,2)] ] ; 
        plot(node_new(1),node_new(2),'g>')
    end
    
    if rem(ij,nnx) == 0 
        node_new = node(ij,:)+[dX/2 0] ;
        node_y = [node_y ; node_new ] ;
        plot(node_new(1),node_new(2),'g>')
    else
        node_new = node(ij,:)+[dX/2 0] ;
        node_y = [node_y ; node_new ] ;
        plot(node_new(1),node_new(2),'b>')
    end
%     pause
end
% plot(node_y(:,1),node_y(:,2),'b>')


