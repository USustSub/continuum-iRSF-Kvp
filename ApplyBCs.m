
    K0 = K; 
    for  ui = 1 : length(rightNodes) 
        cur_node1 = rightNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

    for  ui = 1 : length(botNodes) 
        cur_node1 = botNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

    for  ui = 1 : length(topNodes) 
        cur_node1 = topNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
%         [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

    for  ui = 1 : length(leftNodes) 
        cur_node1 = leftNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
%         [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

%     for  ui = 1 
%         cur_node1 = topNodes ( ui )  ;
%         xx = node(cur_node1,1)/max(node(:,1)); 
%         yy = node(cur_node1,2)/max(node(:,2)); 
%         [K,Residual] = boundary_1point(K,Residual,2*cur_node1,0);
%         [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,0);  
%     end