function [conn_nodes ,node] = find_conn_nodes (node , element ) 

% This function finds all the elements connected to each node

    conn_nodes = zeros(length(node),1) ;  
    cc = zeros(size(node,1),1) ; 
for jj = 1 : size(element,1)
    sctr = element(jj,:) ; 
    if ( sctr(1)==sctr(2)) 
        continue
    end
    for kk = 1 : length(sctr) 
      cc(sctr(kk)) = cc(sctr(kk))+1 ; 
%     sctr(jj,kk) = element(jj,:) ; 
    conn_nodes(sctr(kk),cc(sctr(kk))+1) = jj; 
    conn_nodes(sctr(kk),1) = cc(sctr(kk)); 

    end
end