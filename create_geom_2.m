clear all ; 
for nnx = [80 160 320 380 640]
    nnx
clearvars -except nnx 
close all ;
Ex = 'Preuss' ;
pstep1 = 100 ; % postprocess step 1
pstep2 = 100 ; % postprocess step 2
ifprofile = 0 ; % for debugging
rewrite_mesh = 0 ;

Input_new ; % input parameters and grid generation
Create_geom ;

end