    
function [in,ic,local,o] = find_which_elem ( xx , conn_elems , node , element ) 
    for o = 1 : length(conn_elems)
        ic = conn_elems(o);
        xv = node(element(ic,:),1); xv = [xv;xv(1)];
        yv = node(element(ic,:),2); yv = [yv;yv(1)];
        in = inpolygon(xx(1),xx(2),xv,yv);
        if in == 1
            break
        end
    end
    
        sctr = element(ic,:) ;
        local = element_coordinate ( xx , sctr , node );
