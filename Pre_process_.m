mm_ = 0 ; dY = node(nnx+1,2)-node(1,2);
% dY = dY/2; deltaY = deltaY/2 ; 
DD = [] ;    DD_BC = [] ;           cc_ij = 1 ;
B_x_b = sparse(1*numnode,2*numnode);B_y_b = sparse(1*numnode,2*numnode);
B_xy_b = sparse(1*numnode,2*numnode);B_x_r = sparse(1*numnode,2*numnode);
B_y_r = sparse(1*numnode,2*numnode);B_xy_r = sparse(1*numnode,2*numnode);
B_x_t = sparse(1*numnode,2*numnode);B_y_t = sparse(1*numnode,2*numnode);
B_xy_t = sparse(1*numnode,2*numnode);B_x_l = sparse(1*numnode,2*numnode);
B_y_l = sparse(1*numnode,2*numnode);B_xy_l = sparse(1*numnode,2*numnode);
Dmat_t1_l = zeros(1*numnode,2);Dmat_t1_r = zeros(1*numnode,2);
Dmat_t1_b = zeros(1*numnode,2);Dmat_t1_t = zeros(1*numnode,2);
D11 = zeros(4*numnode,1);D12 = zeros(4*numnode,1);D13 = zeros(4*numnode,1);D14 = zeros(4*numnode,1);
D21 = zeros(4*numnode,1);D22 = zeros(4*numnode,1);D23 = zeros(4*numnode,1);D24 = zeros(4*numnode,1);
D31 = zeros(4*numnode,1);D32 = zeros(4*numnode,1);D33 = zeros(4*numnode,1);D34 = zeros(4*numnode,1);
D41 = zeros(4*numnode,1);D42 = zeros(4*numnode,1);D43 = zeros(4*numnode,1);D44 = zeros(4*numnode,1);
maxS = 0; 
SS = zeros(4*numnode,2) ; 
for ij = 1 : numnode
    ij/numnode
    node_c = node(ij,:) ;

    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;

    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;
    if (abs(node_c(1)-0)>0.0001 && abs(node_c(1)-L)>0.0001  && abs(node_c(2)-D)>0.0001 && abs(node_c(2)-0)>0.0001 )
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
            SS(ij+(ii-1)*numnode,:) = xx ; 

        end
    end
end
SS0 = SS; 
SS(SS(:,1)==0,:)= [];
%{
for ij = 1 : size(node,1)
    ij/size(node,1)
    II = repmat(SS(1:1000,1),[1,length(node)]) ;
end
% dif = node - [ones(numnode,1)*x(1,1) ones(numnode,1)*x(1,2)];
% r = sqrt(sum(abs(dif).^2,2))';
% index = find(r - di <= 0.00001*max(max(node)));
mm_ 

%%
C = (1:10)
M = repmat( C , [10,1] ) ; 
M
%}
%%
id1 = 1 ; 
Bx =  zeros(4*numnode,maxS) ;  
By = zeros(4*numnode,maxS) ; 
En = zeros(4*numnode,2*maxS) ; 

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

            if strcmp(Ex,'Slip') == 1 || strcmp(Ex,'Herr') == 1 
                if norm(abs(xx(2)-D/2))<=dY/2+0.01*dY % && xx(1) > 0.2 && xx(1)<0.8
%                     if norm(abs(xx(1))) < L-dY && norm(abs(xx(1))) > dY
                    DD = [ DD ; xx ij+(ii-1)*numnode] ; 
                        C0 = 0.00 ; 
                        HH(ij,ii*6-0) = C0 ; 
%                     end
                end
            end
            if strcmp(Ex,'Preuss') == 1 
%                 if norm(abs(xx(2)-D/2))<=1.01*dY/2 && norm(abs(xx(1)-L/2))<=dY/2*2+0.00*dY% && xx(1) > 0.2 $
                    rx = xx(1)-L/2;
                    ry = xx(2)-D/2;
                    CL = 1000 ; CL2 = CL/5; 
                    ecc = rx^2/CL/CL + ry^2/CL2/CL2;
                if ecc <= 1
                    DD = [ DD ; xx ij+(ii-1)*numnode] ;
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
                Dmat_t1_l(ij,1:2) = [Gve K_];
            elseif ii == 2
                B_x_r(ij,en) = B_x_r(ij,en) + B(1,:);
                B_y_r(ij,en) = B_y_r(ij,en) + B(2,:);
                B_xy_r(ij,en) = B_xy_r(ij,en) + B(3,:);
                Dmat_t1_r(ij,1:2) = [Gve K_];

            elseif ii == 3
                B_x_b(ij,en) = B_x_b(ij,en) + B(1,:);
                B_y_b(ij,en) = B_y_b(ij,en) + B(2,:);
                B_xy_b(ij,en) = B_xy_b(ij,en) + B(3,:);
                Dmat_t1_b(ij,1:2) = [Gve K_];

            elseif ii == 4
                B_x_t(ij,en) = B_x_t(ij,en) + B(1,:);
                B_y_t(ij,en) = B_y_t(ij,en) + B(2,:);
                B_xy_t(ij,en) = B_xy_t(ij,en) + B(3,:);
                Dmat_t1_t(ij,1:2) = [Gve K_];
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
            
            cc_ij = cc_ij + 1 ; 
        end

    end
            id1 = id1 + 4 ; 
end

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
            if strcmp(Ex,'Preuss')
                P0_ = 1e6*15.2*0.33+[HH(:,1*6-1);HH(:,2*6-1);HH(:,3*6-1);HH(:,4*6-1)];
                HH(:,1*6-1) = 1e6*15.2*0.33 ; 
                HH(:,2*6-1) = 1e6*15.2*0.33 ; 
                HH(:,3*6-1) = 1e6*15.2*0.33 ; 
                HH(:,4*6-1) = 1e6*15.2*0.33 ; 
            end
Gve_ = [ Dmat_t1_l(:,1) ; Dmat_t1_r(:,1); Dmat_t1_b(:,1); Dmat_t1_t(:,1)] ;
K__ = [ Dmat_t1_l(:,2) ; Dmat_t1_r(:,2); Dmat_t1_b(:,2); Dmat_t1_t(:,2)] ;

% shape function info for slip calculations...
%{
xmin_DD = min(DD(:,2))-dY;
xmax_DD = max(DD(:,2))+dY;
maxS_DDt = 0  ; maxS_DDb = 0 ; 
SS_t = [] ; SS_b = [] ; 
for ij = 1 : size(DD,1)
    xij_t = DD(ij,1:2) ; xij_t(2) = xmax_DD;
    xij_b = DD(ij,1:2) ; xij_b(2) = xmin_DD;
    [~, ~,en] = get_data ( xij_t , node , di , form ) ;
    maxS_DDt = max(maxS_DDt,length(en)/2) ; 
    [~, ~,en] = get_data ( xij_b , node , di , form ) ;
    maxS_DDb = max(maxS_DDb,length(en)/2) ; 
    SS_t = [SS_t ; xij_t ]  ;
    SS_b = [SS_b ; xij_b ]  ;
end
phi_DD_t = zeros(size(DD,1),maxS_DDt);
en_DD_t = ones(size(DD,1),maxS_DDt*2);
phi_DD_b = zeros(size(DD,1),maxS_DDb);
en_DD_b = ones(size(DD,1),maxS_DDb*2);
for ij = 1 : size(DD,1)
    xij_t = DD(ij,1:2) ; xij_t(2) = xmax_DD;
    xij_b = DD(ij,1:2) ; xij_b(2) = xmin_DD;
    [phi, ~,en] = get_data ( xij_t , node , di , form ) ;
    phi_DD_t(ij,1:length(phi)) = phi ; 
    en_DD_t(ij,1:length(en)) = en ; 
    [phi, ~,en] = get_data ( xij_b , node , di , form ) ;
    phi_DD_b(ij,1:length(phi)) = phi ; 
    en_DD_b(ij,1:length(en)) = en ; 
end
%}
% save Mesh_data_440.mat 

