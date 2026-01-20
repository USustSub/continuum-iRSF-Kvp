%% old mfree 


K = sparse(2*numnode,2*numnode);
F = zeros(2*numnode,1) ; 

% loop over nodes and assemble stiffness matrix
for ij = 1 : size(node,1)
    if rem(ij,10) == 0 
        ij/size(node,1)
    end
    node_c = node(ij,:) ; 
    
    
    Dmat = identify_tangent ( node_c ) ;

% only for nodes inside domain
if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;    

    get_equations ();
    
end

    get_BCs () ;

end

% get initial elastic solution 
dU = K\F;
% dUx_e    = dU(NumUx);
% dUy_e    = dU(NumUyG);
ddU      = 0*dU;       % Initialise correction vector


fac = 10;
node_deformed = node ; 
node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);

figure
hold on
plot(node(:,1),node(:,2),'bsq')
plot(node_deformed(:,1),node_deformed(:,2),'rsq')
axis equal 

%%

% loop over nodes
    Stress = zeros(numnode,3); 
    Disp = zeros(numnode,2); 
for ij = 1 : size(node,1)
    node_c = node(ij,:) ;
    
    Dmat = identify_tangent ( node_c ) ;
    
    [phi,B,en] = get_data ( node_c , node , di , form ) ;
    Stress(ij,:) = Dmat * B * dU(en) ;
     
    Disp(ij,:) = [phi * dU(en(1:2:end)) phi * dU(en(2:2:end))] ;

end



% tri = delaunay(node(:,1),node(:,2));
% tri = tricheck(node(:,1:2),tri) ; 
tri = element ; 

VTKPostProcess(node,tri,1,'Quad4','stress',Stress,Disp)
!paraview stress.vtu&
 