function [HH,mm_,dY,DD,DD_BC,cc_ij,B_x_b,B_y_b,B_xy_b,B_xy_r,B_y_r,B_x_r,...
B_x_t,B_y_t,B_xy_t,B_x_l,B_y_l,B_xy_l,...
Dmat_t1_l,Dmat_t1_r,Dmat_t1_b,Dmat_t1_t,...
D11,D12,D13,D14,D21,D22,D23,D24,...
D31,D32,D33,D34,D41,D42,D43,D44,...
maxS,SS,SS0,id1,Bx,By,En,Shape,Gve_,K__,dcX , dcY ,dcX_L,dcX_R,dcY_B,dcY_T] = Pre_process_new_efficient_v3 ()
% Get all base workspace variables into this function.
% this version does not use MLS shape functions
T = evalin('base','whos');
for ii = 1:length(T)
   C_ =  evalin('base',[T(ii).name ';']);
   eval([T(ii).name,'=C_;']);
end
clear T C_ ii


mm_ = 0 ; dY = node(nnx+1,2)-node(1,2);

DD = [] ;    DD_BC = [] ;           cc_ij = 1 ;
B_x_b = sparse(1*numnode,2*numnode);B_y_b = sparse(1*numnode,2*numnode);
B_xy_b = sparse(1*numnode,2*numnode);B_x_r = sparse(1*numnode,2*numnode);
B_y_r = sparse(1*numnode,2*numnode);B_xy_r = sparse(1*numnode,2*numnode);
B_x_t = sparse(1*numnode,2*numnode);B_y_t = sparse(1*numnode,2*numnode);
B_xy_t = sparse(1*numnode,2*numnode);B_x_l = sparse(1*numnode,2*numnode);
B_y_l = sparse(1*numnode,2*numnode);B_xy_l = sparse(1*numnode,2*numnode);
SS = zeros(4*numnode,2) ; 
% figure
% hold on
% axis equal
% plot(node(:,1),node(:,2),'bsq')
% disp('evaluating nodes')

dcX = deltaX/2;
dcY = deltaY/2;
dcX = fi_/2 ; 
dcY = fi_/2 ; 
dcX_L = fi_2(:,1)/2;
dcX_R = fi_2(:,2)/2;
dcY_B = fi_2(:,3)/2;
dcY_T = fi_2(:,4)/2;

% dcX = L/1000;  
% dcY = D/1000;  
dcX2 = L/1000 ;
dcY2 = D/1000 ;
dcX2 = deltaX/1000;
dcY2 = deltaY/1000;
    
maxS = 0 ; 
for ig =1:size(conn_nodes,1)
    maxS = max(maxS,length(unique(element(conn_nodes(ig,2:conn_nodes(ig,1)+1),:)))); 
end

id1 = 1 ; 
Bx =  zeros(4*numnode,maxS) ;  
By = zeros(4*numnode,maxS) ; 
En = zeros(4*numnode,2*maxS) ; 
B__l = [] ; B__r = [] ; B__t = [] ; B__b = [] ;

if size(element,2) == 4
    elemType = 'Q4' ; 
else
    elemType = 'T3' ; 
end
for ij = 1 : numnode
%     if rem(ij,round(numnode/10))==0
        disp(['ij -- > ' num2str(ij/numnode*100) ' %'])
%     end
    node_c = node(ij,:) ;

    x_l = node_c - [dcX_L(ij) 0 ] ;
    x_r = node_c + [dcX_R(ij) 0 ] ;

    x_b = node_c - [0 dcY_B(ij) ] ;
    x_t = node_c + [0 dcY_T(ij) ] ;

    sc = length(unique(element(conn_nodes(ij,2:conn_nodes(ij,1)+1),:)));

    if (abs(node_c(1)-0)>0.0001 && abs(node_c(1)-L)>0.0001  && abs(node_c(2)-D)>0.0001 && abs(node_c(2)-0)>0.0001 )
        
            
    conn_elems = conn_nodes(ij,2:conn_nodes(ij,1)+1);
    nn = element(conn_elems,:) ;
    [oe,oe_i,oe_j]  = unique(nn(:)) ;
    oe_2 = unique(oe);
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
%             [index] = define_support(node,xx,di);
            get_index_0 ; 
%             maxS = max(length(index),maxS);
            SS(ij+(ii-1)*numnode,:) = xx ; 
%             plot(node(index,1),node(index,2),'rsq')
%             plot(xx(1),xx(2),'rsq')
            
%     find parent element:
        [in,ic,local,o] = find_which_elem ( xx , conn_elems , node , element ) ;
        xv = node(element(ic,:),1); xv = [xv;xv(1)];
        yv = node(element(ic,:),2); yv = [yv;yv(1)];
%         plot(xv,yv,xx(in,1),xx(in,2),'.r');
%         pause

        xx_l = xx - [dcX2 0 ] ;
        [~,ic_l,local_l,~] = find_which_elem ( xx_l , conn_elems , node , element ) ;
        xx_r = xx + [dcX2 0 ] ;
        [~,ic_r,local_r,~] = find_which_elem ( xx_r , conn_elems , node , element ) ;
        xx_b = xx - [0 dcY2 ] ;
        [~,ic_b,local_b,~] = find_which_elem ( xx_b , conn_elems , node , element ) ;
        xx_t = xx + [0 dcY2 ] ;
        [~,ic_t,local_t,~] = find_which_elem ( xx_t , conn_elems , node , element ) ;

%         get_index ;             
        index2 = oe' ;
        en = zeros(1,2*size(index2,2));  
        for m = 1 : size(index2,2)      
            en(2*m-1) = 2*index2(m)-1;
            en(2*m  ) = 2*index2(m)  ;
        end    

% find B
        [Nv_l,~]=lagrange_basis(elemType,local_l,1) ;
        [Nv_r,~]=lagrange_basis(elemType,local_r,1) ;
        [Nv_t,~]=lagrange_basis(elemType,local_t,1) ;
        [Nv_b,~]=lagrange_basis(elemType,local_b,1) ;
        dphidx = zeros(length(node),1) ; 
        dphidx(element(ic_r,:)) = dphidx(element(ic_r,:)) + Nv_r/(dcX2*2) ;
        dphidx(element(ic_l,:)) = dphidx(element(ic_l,:)) - Nv_l/(dcX2*2) ;
        dphidx = dphidx(index2);
        
        dphidy = zeros(length(node),1) ; 
        dphidy(element(ic_t,:)) = dphidy(element(ic_t,:)) + Nv_t/(dcY2*2) ;
        dphidy(element(ic_b,:)) = dphidy(element(ic_b,:)) - Nv_b/(dcY2*2) ;
        dphidy = dphidy(index2);
        
        B = zeros(3,sc*2) ; 
        B(1,1:2:end) = dphidx ;
        B(2,2:2:end) = dphidy ;
        B(3,1:2:end) = dphidy ; 
        B(3,2:2:end) = dphidx ;
    


        dNdx = zeros(1,maxS) ; dNdx(:) = 0*-10000.1 ; 
        dNdy = zeros(1,maxS) ; dNdy(:) = 0*-10000.1 ;  
        en_ = zeros(1,2*maxS) ; en_(:) = -10000.1 ;  
        dNdx(1:size(B,2)/2) = B(1,1:2:end);
        dNdy(1:size(B,2)/2) = B(2,2:2:end);
        en_(1:size(en,2)) = en;
    
        Bx(ij+(ii-1)*numnode,:) = [ dNdx ] ;
        By(ij+(ii-1)*numnode,:) = [ dNdy ] ;
        En(ij+(ii-1)*numnode,:) = [ en_ ] ;
            
            

        if ii == 1
            B_x_l(ij,en) = B_x_l(ij,en) + B(1,:);
            B_y_l(ij,en) = B_y_l(ij,en) + B(2,:);
            B_xy_l(ij,en) = B_xy_l(ij,en) + B(3,:);
        elseif ii == 2
            B_x_r(ij,en) = B_x_r(ij,en) + B(1,:);
            B_y_r(ij,en) = B_y_r(ij,en) + B(2,:);
            B_xy_r(ij,en) = B_xy_r(ij,en) + B(3,:);
        elseif ii == 3
            B_x_b(ij,en) = B_x_b(ij,en) + B(1,:);
            B_y_b(ij,en) = B_y_b(ij,en) + B(2,:);
            B_xy_b(ij,en) = B_xy_b(ij,en) + B(3,:);
        elseif ii == 4
            B_x_t(ij,en) = B_x_t(ij,en) + B(1,:);
            B_y_t(ij,en) = B_y_t(ij,en) + B(2,:);
            B_xy_t(ij,en) = B_xy_t(ij,en) + B(3,:);
        end
            
            Shape{cc_ij}{1} = B;
            Shape{cc_ij}{2} = en ;
            
            cc_ij = cc_ij + 1 ; 
        

        end
    end
end
SS0 = SS; 
SS(SS(:,1)==0,:)= [];



%% get material properties
HH = zeros(size(node,1),(4+1+1)*4); % history parameter for L(1), B(2), R(3), T(4)
HH(:,[6 12 18 24]) = coh0;

Dmat_t1_l = zeros(1*numnode,2);Dmat_t1_r = zeros(1*numnode,2);
Dmat_t1_b = zeros(1*numnode,2);Dmat_t1_t = zeros(1*numnode,2);
D11 = zeros(4*numnode,1);D12 = zeros(4*numnode,1);D13 = zeros(4*numnode,1);D14 = zeros(4*numnode,1);
D21 = zeros(4*numnode,1);D22 = zeros(4*numnode,1);D23 = zeros(4*numnode,1);D24 = zeros(4*numnode,1);
D31 = zeros(4*numnode,1);D32 = zeros(4*numnode,1);D33 = zeros(4*numnode,1);D34 = zeros(4*numnode,1);
D41 = zeros(4*numnode,1);D42 = zeros(4*numnode,1);D43 = zeros(4*numnode,1);D44 = zeros(4*numnode,1);

disp('evaluating material properties')
for ij = 1 : numnode % loop over nodes 
    if rem(ij,numnode/10)==0
    disp(['ij -- > ' num2str(ij/numnode*100) ' %'])
    end
    node_c = node(ij,:) ;

    x_l = node_c - [dcX_L(ij) 0 ] ;
    x_r = node_c + [dcX_R(ij) 0 ] ;

    x_b = node_c - [0 dcY_B(ij) ] ;
    x_t = node_c + [0 dcY_T(ij) ] ;

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

            [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( xx , mat , dynamic , node) ;
            

            if ii == 1
                Dmat_t1_l(ij,1:2) = [Gve K_];
            elseif ii == 2
                Dmat_t1_r(ij,1:2) = [Gve K_];

            elseif ii == 3
                Dmat_t1_b(ij,1:2) = [Gve K_];

            elseif ii == 4
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
        end

    end
            id1 = id1 + 4 ; 
end


Gve_ = [ Dmat_t1_l(:,1) ; Dmat_t1_r(:,1); Dmat_t1_b(:,1); Dmat_t1_t(:,1)] ;
K__ = [ Dmat_t1_l(:,2) ; Dmat_t1_r(:,2); Dmat_t1_b(:,2); Dmat_t1_t(:,2)] ;




