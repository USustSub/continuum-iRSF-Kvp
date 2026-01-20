theta_ = 0 ; 

% ! rm -r in/
% ! mkdir in/

! rm -r stress
! mkdir stress
! rm -r out
! mkdir out
! mkdir out/paraview
! mkdir out/slip
! mkdir out/slip_center
! mkdir out/slip_rate
! mkdir out/theta
! mkdir out/mu
! mkdir out/deformed
! mkdir out/dlambda

! mkdir out/paraview/theta
! mkdir out/paraview/mu
! mkdir out/paraview/sxx
! mkdir out/paraview/syy
! mkdir out/paraview/szz
! mkdir out/paraview/sxy
! mkdir out/paraview/dggg
! mkdir out/paraview/Eii
! mkdir out/paraview/ax
! mkdir out/paraview/vx
! mkdir out/paraview/P
! mkdir out/paraview/J
! mkdir out/history

! for dir in out/*/; do mkdir -- "$dir/data"; done

warning('off','MATLAB:nearlySingularMatrix');
addpath ('/home/mohsen/Downloads/XFEMMatlabcode')
addpath ('/home/mohsen/Desktop/MHSN/MatlabTools/p1top2')


inc       = 0; % number of increment
ninc          = 20000000 ;     % Number of increments
gitmax        = 2500;     % Max. number of global iterations
tol_plast     = 1e-10;  % Tolerance of global iterations


if strcmp(Ex,'Duretz')
    input_duretz ;
elseif strcmp(Ex,'Slip')
    input_slip ;
elseif strcmp(Ex,'Herr') || strcmp(Ex,'Junction')
    input_H ; 
elseif strcmp(Ex,'Preuss')
    input_P ;
elseif strcmp(Ex,'Preuss2')
    input_P ;
end
