%% 2D MLS based FDM 
% By Mohsen Goudarzi ( m.goudarzi@uu.nl ) 
% ----------------------------------------------------------
clear all; close all ; 
clc
tic;           % help us to see the time required for each step

% Dimension of the domain (it is simply a rectangular region L x W)
L = 4 ;
D = 4 ;

% Material properties
E  = 100 ;
nu = 0.3 ;


% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 2.0 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function


% Node density defined by number of nodes along two directions
nnx = 20 ;
nny = 20 ;

% node generation, node is of size (nnode x 2)
node = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);
numnode = size(node,1);

deltaX = L/(nnx-1);
deltaY = D/(nny-1);
delta  = max(deltaX,deltaY);
di     = ones(1,numnode)*dmax*delta ;



% Meshing, Q4 elements
inc_u = 1;
inc_v = nnx;
node_pattern = [ 1 2 nnx+2 nnx+1 ];
element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);

%

figure
hold on
plot(node(:,1),node(:,2),'rsq')
for icel = 1 : size(element,1)
    sc = element(icel,:) ; 
    node_c = mean(node(sc,:));
    
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    y_b = node_c - [0 deltaX/2 ] ;
    y_t = node_c + [0 deltaX/2 ] ;    

end


%% add crack 
a = 2 ; 
xCr = [ 0 D/2 ; a D/2 ] ;
xTip = [a D/2];
seg  = xCr(2,:) - xCr(1,:);   % tip segment


x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);
t   = 1/norm(seg)*seg;
for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip)*t';  % tangent LS
    dt(i)   = norm([x y]-xTip); % distances from node I to tip
end

set1 = find(abs(ls(:,1)) - di' < 0);
set2 = find(ls(:,2) < 0);
set3 = find(dt - di < 0);        % tip nodes
set4 = intersect(set1,set2);     % split nodes

% some nodes belong to both sets, remove tip
split_nodes = setdiff(set4,set3);

plot(node(split_nodes,1),node(split_nodes,2),'bsq')

pos = zeros(numnode,1);
nsnode = 0 ;
for i = 1 : numnode
    [snode] = ismember(i,split_nodes);
    if (snode ~= 0)
        pos(i) = (numnode + nsnode*1 ) + 1 ;
        nsnode = nsnode + 1 ;
    end
end
tip_nodes = [] ;
 
%%
Consti;

total_num_node = numnode + length(split_nodes);
K = sparse(2*total_num_node,2*total_num_node);
F = zeros(2*total_num_node,1);

% loop over nodes
for ij = 1 : size(node,1)
    node_c = node(ij,:) ; 
    
if (node_c(1) > 0 && node_c(1) < L && node_c(2) < L && node_c(2) > 0 )
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;    

% x momentum equation - left
    [~,B,en] = get_data_enr ( x_l , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;

% x momentum equation - right
    [~,B,en] = get_data_enr ( x_r , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% x momentum equation - bot
    [~,B,en] = get_data_enr ( x_b , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;

% x momentum equation - top
    [~,B,en] = get_data_enr ( x_t , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% y momentum equation - left
    [~,B,en] = get_data_enr ( x_l , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - right
    [~,B,en] = get_data_enr ( x_r , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;

% y momentum equation - bot
    [~,B,en] = get_data_enr ( x_b , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - top
    [~,B,en] = get_data_enr ( x_t , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;

end


if node_c(2) == 0 % bot edge
       plot(node_c(1),node_c(2),'bsq')
% only drichlet boundary conditions
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
elseif node_c(2) == L % top edge
       plot(node_c(1),node_c(2),'ksq')
% drichlet (ux=0.05) and neumann (t2=0) boundary conditions
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)/L*0.05 ;
        
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
end


% Neuman boundary condition
if node_c(1) == 0 && node_c(2) >0 && node_c(2) < L  % left edge
% neumann (t1=t2=0) boundary conditions
       plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(1,:)*B; %+ Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
elseif node_c(1) == L && node_c(2) >0 && node_c(2) < L  % left edge
% neumann (t1=t2=0) boundary conditions
       plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(1,:)*B; % ; + Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
end

 
end


% loop over enrichment nodes
for ij_ = 1 : length(split_nodes)
    node_c = node(split_nodes(ij_),:) ; 
    ij = pos(split_nodes(ij_)) ; 
    
if (node_c(1) > 0 && node_c(1) < L && node_c(2) < L && node_c(2) > 0 )
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;    

% x momentum equation - left
    [~,B,en] = get_data_enr ( x_l , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;

% x momentum equation - right
    [~,B,en] = get_data_enr ( x_r , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% x momentum equation - bot
    [~,B,en] = get_data_enr ( x_b , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;

% x momentum equation - top
    [~,B,en] = get_data_enr ( x_t , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% y momentum equation - left
    [~,B,en] = get_data_enr ( x_l , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - right
    [~,B,en] = get_data_enr ( x_r , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;

% y momentum equation - bot
    [~,B,en] = get_data_enr ( x_b , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - top
    [~,B,en] = get_data_enr ( x_t , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;
end


if node_c(2) == 0 % bot edge
       plot(node_c(1),node_c(2),'bsq')
% only drichlet boundary conditions
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
elseif node_c(2) == L % top edge
       plot(node_c(1),node_c(2),'ksq')
% drichlet (ux=0.05) and neumann (t2=0) boundary conditions
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)/L*0.05 ;
        
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
end


% Neuman boundary condition
if node_c(1) == 0 && node_c(2) >0 && node_c(2) < L  % left edge
% neumann (t1=t2=0) boundary conditions
       plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(1,:)*B; %+ Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
elseif node_c(1) == L && node_c(2) >0 && node_c(2) < L  % left edge
% neumann (t1=t2=0) boundary conditions
       plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(1,:)*B; % ; + Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
end

 
end


figure
spy(K)



%% apply traction free boundary condition on the crack surface (todo)
ss = [0 :  0.05 : xCr(1,2)];



return
%%


u = K\F ; 

fac = 1;
node_deformed = node ; 
node_deformed(:,1) = node(:,1) + fac*u(1:2:end);
node_deformed(:,2) = node(:,2) + fac*u(2:2:end);

figure
hold on
plot(node(:,1),node(:,2),'bsq')
plot(node_deformed(:,1),node_deformed(:,2),'rsq')
axis equal 



%% compute stress at nodes


% loop over nodes
    Stress = zeros(numnode,3); 
    Disp = zeros(numnode,2); 
for ij = 1 : size(node,1)
    node_c = node(ij,:) ; 
    [phi,B,en] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) ;
    Stress(ij,:) = Dmat * B * u(en) ;
     
    Disp(ij,:) = [phi * u(en(1:2:end)) phi * u(en(2:2:end))] ;

    
end



tri = delaunay(node(:,1),node(:,2));
% tri = tricheck(node(:,1:2),tri) ; 

VTKPostProcess(node,tri,1,'Tri3','stress',Stress,Disp)
!paraview stress.vtu&







