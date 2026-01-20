function [sctrB,snode,tnode] = assembly(index,split_nodes,tip_nodes,...
                                        pos)

[snode,s_loc] = ismember(index,split_nodes);
[tnode,t_loc] = ismember(index,tip_nodes);

% number of H and tip enriched nodes 
% in index
num_split_node = size(find(snode ~= 0),2);
num_tip_node   = size(find(tnode ~= 0),2);

% Scatter of standard part of B matrix
le       = size(index,2);
sctrStdB = zeros(2*le,1);
sctrStdB(1:2:2*le) = index.*2-1 ; % x displacement
sctrStdB(2:2:2*le) = index.*2   ; % y displacement

% No enriched nodes in index
if (num_split_node == 0) && (num_tip_node == 0)
    sctrB = sctrStdB ;
    % There is at least one enriched node in index
else
    sn = num_split_node;
    tn = num_tip_node  ;
    sctrEnrB = zeros(2*(sn*1+tn*4),1);
    cnt = 0 ;
    for k = 1 : le
        nodeI = index(k) ; % I-th node in index
        if (snode(k) ~= 0)
            cnt = cnt + 1 ;
            sctrEnrB(2*cnt - 1) = 2 * pos(nodeI) - 1;
            sctrEnrB(2*cnt    ) = 2 * pos(nodeI)    ;
        elseif (tnode(k) ~= 0)
            cnt = cnt + 1 ;
            sctrEnrB(2*cnt - 1) = 2 * pos(nodeI) - 1;
            sctrEnrB(2*cnt    ) = 2 * pos(nodeI)    ;

            cnt = cnt + 1 ;
            sctrEnrB(2*cnt - 1) = 2 * (pos(nodeI)+1) - 1;
            sctrEnrB(2*cnt    ) = 2 * (pos(nodeI)+1)    ;

            cnt = cnt + 1 ;
            sctrEnrB(2*cnt - 1) = 2 * (pos(nodeI)+2) - 1;
            sctrEnrB(2*cnt    ) = 2 * (pos(nodeI)+2)    ;

            cnt = cnt + 1 ;
            sctrEnrB(2*cnt - 1) = 2 * (pos(nodeI)+3) - 1;
            sctrEnrB(2*cnt    ) = 2 * (pos(nodeI)+3)    ;
        end
    end
    sctrB = [ sctrStdB;sctrEnrB ];
end