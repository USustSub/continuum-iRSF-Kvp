
function [phi, B , sctrB] = get_data_enr ( node_c , node , di , form , split_nodes,tip_nodes,pos,xCr,xTip) 
    [index] = define_support(node,node_c,di);
    B = zeros(3,2*size(index,2)) ;
    en = zeros(1,2*size(index,2));  
  
    [phi,dphidx,dphidy] = MLS_ShapeFunction(node_c,index,node,di,form);
    for m = 1 : size(index,2)      
        B(1:3,2*m-1:2*m) = [dphidx(m) 0 ; 
                           0 dphidy(m);
                           dphidy(m) dphidx(m)];
        en(2*m-1) = 2*index(m)-1;
        en(2*m  ) = 2*index(m)  ;
    end    

[sctrB,snode,tnode] = assembly(index,split_nodes,tip_nodes,pos);

    %  compute B matrix
   [ B , phi ] = Bmatrix(node_c,index,node,di,form,...
                     snode,tnode,xCr,xTip,0);
