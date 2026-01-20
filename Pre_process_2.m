mm_ = 0 ; dY = node(nnx+1,2)-node(1,2);
DD = [] ;   DD_BC = [] ;      cc_ij = 1 ;
B_x_b = sparse(1*numnode,2*numnode);B_y_b = sparse(1*numnode,2*numnode);
B_xy_b = sparse(1*numnode,2*numnode);B_x_r = sparse(1*numnode,2*numnode);
B_y_r = sparse(1*numnode,2*numnode);B_xy_r = sparse(1*numnode,2*numnode);
B_x_t = sparse(1*numnode,2*numnode);B_y_t = sparse(1*numnode,2*numnode);
B_xy_t = sparse(1*numnode,2*numnode);B_x_l = sparse(1*numnode,2*numnode);
B_y_l = sparse(1*numnode,2*numnode);B_xy_l = sparse(1*numnode,2*numnode);
Dmat_t1_l = sparse(1*numnode,6);Dmat_t1_r = sparse(1*numnode,6);
Dmat_t1_b = sparse(1*numnode,6);Dmat_t1_t = sparse(1*numnode,6);
Dmat_t2_l = sparse(1*numnode,6);Dmat_t2_r = sparse(1*numnode,6);
Dmat_t2_b = sparse(1*numnode,6);Dmat_t2_t = sparse(1*numnode,6);
Dmat_t3_l = sparse(1*numnode,6);Dmat_t3_r = sparse(1*numnode,6);
Dmat_t3_b = sparse(1*numnode,6);Dmat_t3_t = sparse(1*numnode,6);
Dmat_t4_l = sparse(1*numnode,6);Dmat_t4_r = sparse(1*numnode,6);
Dmat_t4_b = sparse(1*numnode,6);Dmat_t4_t = sparse(1*numnode,6);
Dmat_t1_l2  = sparse(1*numnode,4*numnode);
Dmat_t2_l2  = sparse(1*numnode,4*numnode);
Dmat_t3_l2  = sparse(1*numnode,4*numnode);
Dmat_t4_l2  = sparse(1*numnode,4*numnode);
Dmat_t1_r2  = sparse(1*numnode,4*numnode);
Dmat_t2_r2  = sparse(1*numnode,4*numnode);
Dmat_t3_r2  = sparse(1*numnode,4*numnode);
Dmat_t4_r2  = sparse(1*numnode,4*numnode);
Dmat_t1_b2  = sparse(1*numnode,4*numnode);
Dmat_t2_b2  = sparse(1*numnode,4*numnode);
Dmat_t3_b2  = sparse(1*numnode,4*numnode);
Dmat_t4_b2  = sparse(1*numnode,4*numnode);
Dmat_t1_t2  = sparse(1*numnode,4*numnode);
Dmat_t2_t2  = sparse(1*numnode,4*numnode);
Dmat_t3_t2  = sparse(1*numnode,4*numnode);
Dmat_t4_t2  = sparse(1*numnode,4*numnode);
D11 = zeros(4*numnode,1);D12 = zeros(4*numnode,1);D13 = zeros(4*numnode,1);D14 = zeros(4*numnode,1);
D21 = zeros(4*numnode,1);D22 = zeros(4*numnode,1);D23 = zeros(4*numnode,1);D24 = zeros(4*numnode,1);
D31 = zeros(4*numnode,1);D32 = zeros(4*numnode,1);D33 = zeros(4*numnode,1);D34 = zeros(4*numnode,1);
D41 = zeros(4*numnode,1);D42 = zeros(4*numnode,1);D43 = zeros(4*numnode,1);D44 = zeros(4*numnode,1);
maxS = 0; 

for ij = 1 : numnode
    ij/numnode
    node_c = node(ij,:) ;

    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;

    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;
 
    if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
        for ii = 1 : 4
            if ii == 1
                xx = x_l ; 
            elseif ii == 2
                xx = x_r ;
            elseif ii == 3
                xx = x_b ;
            elseif ii == 4
                xx = x_t ;
            end
            mm_ = mm_  + 1 ;  
            [index] = define_support(node,xx,di);
            maxS = max(length(index),maxS);
        end
    end
end
mm_
%%
id1 = 1 ; 
Bx =  sparse(4*numnode,maxS) ;  
By = sparse(4*numnode,maxS) ; 
En = sparse(4*numnode,2*maxS) ; 
SS = zeros(4*numnode,2) ; 

for ij = 1 : numnode % loop over nodes 
    ij/numnode
    node_c = node(ij,:) ;

    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;

    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;

    if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
        for ii = 1 : 4
            if ii == 1
                xx = x_l ; 
            elseif ii == 2
                xx = x_r ;
            elseif ii == 3
                xx = x_b ;
            elseif ii == 4
                xx = x_t ;
            end

            id1 = ij+ii*numnode-numnode ;
            if strcmp(Ex,'Slip') == 1 || strcmp(Ex,'Herr') == 1 
                if norm(abs(xx(2)-max(node(:,2))/2))<dY % && xx(1) > 0.2 && xx(1)<0.8
                    DD = [ DD ; xx id1] ; 
                        C0 = 0.00 ; 
                        HH(ij,ii*6-0) = C0 ; 
                end
            end
            
            [~,B,en] = get_data ( xx , node , di , form ) ;
            dNdx = zeros(1,maxS) ; dNdx(:) = 0*-10000.1 ; 
            dNdy = zeros(1,maxS) ; dNdy(:) = 0*-10000.1 ;  
            en_ = zeros(1,2*maxS) ; en_(:) = -10000.1 ;  
            dNdx(1:size(B,2)/2) = B(1,1:2:end);
            dNdy(1:size(B,2)/2) = B(2,2:2:end);
            en_(1:size(en,2)) = en;
    
            Bx(ij+(ii-1)*numnode,:) = [ dNdx ] ;
            By(ij+(ii-1)*numnode,:) = [ dNdy ] ;
            En(ij+(ii-1)*numnode,:) = [ en_ ] ;
            
            [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( xx , mat , dynamic , node) ;
            if ii == 1
                B_x_l(ij,en) = B_x_l(ij,en) + B(1,:);
                B_y_l(ij,en) = B_y_l(ij,en) + B(2,:);
                B_xy_l(ij,en) = B_xy_l(ij,en) + B(3,:);
                Dmat_t1_l(ij,1:4) = Dmat_t(1,:); Dmat_t1_l(ij,5:6) = [Gve K_];
                Dmat_t2_l(ij,1:4) = Dmat_t(2,:); Dmat_t2_l(ij,5:6) = [Gve K_];
                Dmat_t3_l(ij,1:4) = Dmat_t(3,:); Dmat_t3_l(ij,5:6) = [Gve K_]; 
                Dmat_t4_l(ij,1:4) = Dmat_t(4,:); Dmat_t4_l(ij,5:6) = [Gve K_]; 
                Dmat_t1_l2(ij,id1:id1+3) = Dmat_t(1,1:4) ;
                Dmat_t2_l2(ij,id1:id1+3) = Dmat_t(2,1:4) ;
                Dmat_t3_l2(ij,id1:id1+3) = Dmat_t(3,1:4) ;
                Dmat_t4_l2(ij,id1:id1+3) = Dmat_t(4,1:4) ;

            elseif ii == 2
                B_x_r(ij,en) = B_x_r(ij,en) + B(1,:);
                B_y_r(ij,en) = B_y_r(ij,en) + B(2,:);
                B_xy_r(ij,en) = B_xy_r(ij,en) + B(3,:);
                Dmat_t1_r(ij,1:4) = Dmat_t(1,:); Dmat_t1_r(ij,5:6) = [Gve K_];
                Dmat_t2_r(ij,1:4) = Dmat_t(2,:); Dmat_t2_r(ij,5:6) = [Gve K_];
                Dmat_t3_r(ij,1:4) = Dmat_t(3,:); Dmat_t3_r(ij,5:6) = [Gve K_]; 
                Dmat_t4_r(ij,1:4) = Dmat_t(4,:); Dmat_t4_r(ij,5:6) = [Gve K_]; 
                Dmat_t1_r2(ij,id1:id1+3) = Dmat_t(1,1:4) ;
                Dmat_t2_r2(ij,id1:id1+3) = Dmat_t(2,1:4) ;
                Dmat_t3_r2(ij,id1:id1+3) = Dmat_t(3,1:4) ;
                Dmat_t4_r2(ij,id1:id1+3) = Dmat_t(4,1:4) ;
            elseif ii == 3
                B_x_b(ij,en) = B_x_b(ij,en) + B(1,:);
                B_y_b(ij,en) = B_y_b(ij,en) + B(2,:);
                B_xy_b(ij,en) = B_xy_b(ij,en) + B(3,:);
                Dmat_t1_b(ij,1:4) = Dmat_t(1,:); Dmat_t1_b(ij,5:6) = [Gve K_];
                Dmat_t2_b(ij,1:4) = Dmat_t(2,:); Dmat_t2_b(ij,5:6) = [Gve K_];
                Dmat_t3_b(ij,1:4) = Dmat_t(3,:); Dmat_t3_b(ij,5:6) = [Gve K_]; 
                Dmat_t4_b(ij,1:4) = Dmat_t(4,:); Dmat_t4_b(ij,5:6) = [Gve K_]; 
                Dmat_t1_b2(ij,id1:id1+3) = Dmat_t(1,1:4) ;
                Dmat_t2_b2(ij,id1:id1+3) = Dmat_t(2,1:4) ;
                Dmat_t3_b2(ij,id1:id1+3) = Dmat_t(3,1:4) ;
                Dmat_t4_b2(ij,id1:id1+3) = Dmat_t(4,1:4) ;
            elseif ii == 4
                B_x_t(ij,en) = B_x_t(ij,en) + B(1,:);
                B_y_t(ij,en) = B_y_t(ij,en) + B(2,:);
                B_xy_t(ij,en) = B_xy_t(ij,en) + B(3,:);
                Dmat_t1_t(ij,1:4) = Dmat_t(1,:); Dmat_t1_t(ij,5:6) = [Gve K_];
                Dmat_t2_t(ij,1:4) = Dmat_t(2,:); Dmat_t2_t(ij,5:6) = [Gve K_];
                Dmat_t3_t(ij,1:4) = Dmat_t(3,:); Dmat_t3_t(ij,5:6) = [Gve K_]; 
                Dmat_t4_t(ij,1:4) = Dmat_t(4,:); Dmat_t4_t(ij,5:6) = [Gve K_]; 
                Dmat_t1_t2(ij,id1:id1+3) = Dmat_t(1,1:4) ;
                Dmat_t2_t2(ij,id1:id1+3) = Dmat_t(2,1:4) ;
                Dmat_t3_t2(ij,id1:id1+3) = Dmat_t(3,1:4) ;
                Dmat_t4_t2(ij,id1:id1+3) = Dmat_t(4,1:4) ;
            end
            
            
            D11(ij+(ii-1)*numnode)=Dmat_t(1,1);
            D12(ij+(ii-1)*numnode)=Dmat_t(1,2);
            D13(ij+(ii-1)*numnode)=Dmat_t(1,3);
            D14(ij+(ii-1)*numnode)=Dmat_t(1,4);
            D21(ij+(ii-1)*numnode)=Dmat_t(2,1);
            D22(ij+(ii-1)*numnode)=Dmat_t(2,2);
            D23(ij+(ii-1)*numnode)=Dmat_t(2,3);
            D24(ij+(ii-1)*numnode)=Dmat_t(2,4);
            D31(ij+(ii-1)*numnode)=Dmat_t(3,1);
            D32(ij+(ii-1)*numnode)=Dmat_t(3,2);
            D33(ij+(ii-1)*numnode)=Dmat_t(3,3);
            D34(ij+(ii-1)*numnode)=Dmat_t(3,4);
            D41(ij+(ii-1)*numnode)=Dmat_t(4,1);
            D42(ij+(ii-1)*numnode)=Dmat_t(4,2);
            D43(ij+(ii-1)*numnode)=Dmat_t(4,3);
            D44(ij+(ii-1)*numnode)=Dmat_t(4,4);
            
            Shape{cc_ij}{1} = B;
            Shape{cc_ij}{2} = en ;
            
            SS(ij+(ii-1)*numnode,:) = xx ; 
            cc_ij = cc_ij + 1 ; 
        end

    end
            id1 = id1 + 4 ; 
end
SS0 = SS; 
SS(SS(:,1)==0,:)= [];

% initial C
            C0 = [HH(:,1*6-0);HH(:,2*6-0);HH(:,3*6-0);HH(:,4*6-0)];
           
% initial P
            if strcmp(Ex,'Herr')
                P0_ = 5e6+[HH(:,1*6-1);HH(:,2*6-1);HH(:,3*6-1);HH(:,4*6-1)];
                HH(:,1*6-1) = 5e6 ; 
                HH(:,2*6-1) = 5e6 ; 
                HH(:,3*6-1) = 5e6 ; 
                HH(:,4*6-1) = 5e6 ; 
            end
Gve_ = [ Dmat_t1_l(:,5) ; Dmat_t1_r(:,5); Dmat_t1_b(:,5); Dmat_t1_t(:,5)] ;
K__ = [ Dmat_t1_l(:,6) ; Dmat_t1_r(:,6); Dmat_t1_b(:,6); Dmat_t1_t(:,6)] ;

% save Mesh_data_440.mat 