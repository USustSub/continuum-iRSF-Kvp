    
%% Solve system   
	numnode = length(node) ; 
	dispNodes1 = botNodes;
	dispNodes2 = topNodes; 
	dispNodes3 = [botNodes;topNodes;rightNodes;leftNodes] ; 
    
    anode = 1:2*(numnode);    % index vector of all DOF's
    dispNodes = [dispNodes1;dispNodes2;dispNodes3];
    presDofs = [dispNodes1*2-1 ;dispNodes2*2-1 ; dispNodes3*2 ];

%   total prescribed disp vector 
    up = 0*anode'; 
    up (dispNodes1*2-1) = ap2 ; 
    up (dispNodes2*2-1) = ap1 ; 
%  prescribed values 
    ap = zeros(length(presDofs),1) ; 
    ap(1:length(dispNodes1))=ap2 ; 
    ap(length(dispNodes1)+1:length(dispNodes1)+length(dispNodes2)) = ap1 ;
    A = K + a0*Mass ;
    
%     Converged_last_step = Mass*a0*dU0 + Mass*a2*vu0 + Mass*a4*au0 ;

%     r = [Residual+dynamic*Converged_last_step] ;
    r = Residual ; 
    [A,r] = Essential_FEM (A,r,dispNodes,presDofs,anode,plast_it,ap);
    dd = A\r;
    


%     du = dd(1:DOF) ;
    du = 0*up ;
    du(setdiff(anode,presDofs))=dd ; 

% update solution
    dU = dU - du ;
 
% replace actual values for essential boundaries
    dU (dispNodes1*2-1) = dU0 (dispNodes1*2-1) + ap2 ; 
    dU (dispNodes2*2-1) = dU0 (dispNodes2*2-1) + ap1 ; 

    nR = norm(du);
%    if rem(itime,500)==0
    disp ([ 'normm (' num2str(plast_it) ') = ' num2str(nR) ' at  itime = ' num2str(inc) ' maxF --> ' num2str(maxF)  ]); 
