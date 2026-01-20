lambda_matrix = nu*E/( (1+nu)*(1-2*nu) );
mu_matrix     = E/( 2*(1+nu) ); 
            
Dmat = zeros(3);
Dmat(1,1) = lambda_matrix + 2*mu_matrix; 
Dmat(2,2) = lambda_matrix + 2*mu_matrix; 
Dmat(3,3) = mu_matrix; 
Dmat(1,2) = lambda_matrix;
Dmat(2,1) = lambda_matrix; 