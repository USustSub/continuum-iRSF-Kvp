function [phi, B,en] = get_data_BC ( x_l , node , di , form ) 
    [index] = define_support(node,x_l,di);
    B = zeros(3,2*size(index,2)) ;
    en = zeros(1,2*size(index,2));  
  
    [phi,dphidx,dphidy] = MLS_ShapeFunction(x_l,index,node,di,form);
    phi(abs(phi)<0.5)=0;
    phi(abs(phi)>0.5)=1;
    dphidx(abs(dphidx)<1e-4)=0;
    dphidy(abs(dphidy)<1e-4)=0;
    dphidx = sign(dphidx);
    dphidy = sign(dphidy);
% %     dphidx = dphidx/max(abs(dphidx));
% %     dphidx_t = 0*dphidx ; 
% %     dphidx_t(dphidx<-0.1)=-1;
% %     dphidx_t(dphidx>0.1)=1;
% %     dphidx = dphidx_t; 
% %     
% %     dphidy = dphidy/max(abs(dphidy));
% %     dphidy_t = 0*dphidy ; 
% %     dphidy_t(dphidy<-0.1)=-1;
% %     dphidy_t(dphidy>0.1)=1;
% %     dphidy = dphidy_t; 
% %     
    
    if any(isnan(phi))
        phi
    end
    
    for m = 1 : size(index,2)      
        B(1:3,2*m-1:2*m) = [dphidx(m) 0 ; 
                           0 dphidy(m);
                           dphidy(m) dphidx(m)];
        en(2*m-1) = 2*index(m)-1;
        en(2*m  ) = 2*index(m)  ;
    end    

