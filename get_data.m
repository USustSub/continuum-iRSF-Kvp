function [phi, B,en] = get_data ( x_l , node , di , form ) 
    [index] = define_support(node,x_l,di);
    B = zeros(3,2*size(index,2)) ;
    en = zeros(1,2*size(index,2));  
  
    [phi,dphidx,dphidy] = MLS_ShapeFunction(x_l,index,node,di,form);
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

