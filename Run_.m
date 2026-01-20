% run for different diameters and mesh densities

clear all 

! rm -r out_results/
! mkdir out_results/

% problem = 4000; 
% eval(['!sed -i ',char(39),'32s/.*/problem = ',num2str(problem),...
%     ' /',char(39),' Main.m']);


% eval(['!sed -i ',char(39),'3s/.*/L = ',num2str(L,10),...
%      '; /',char(39),' in/3SolidFib.geo']);
% for ii__ = [11 31 61]
% for ii__ = [ 2.5e5 2.5e15 2.5e16 2.5e17 2.5e18 ]

% for df__ = [0.01 0.05 0.1 0.5 1 ]
L = 18750 ; 
for df__ = [ L/100  ]
    
% for ii__ = 112*[ 16 2 4 8 ]% 840 1280]
for ii__ = 60*[ 4 8 16]% 840 1280]
% for nf__ = [ 200 500 1000 10000 20000 ]
for nf__ = 1
    

nF = [ 0 2.5e8 2.5e13 2.5e16 0 2.5e17 2.5e18 ] ; 

ee__ = nF(nf__) ; 
% 
% tend = fliplr([2000 1000 500 200 100]);    
% tend(dd__)
    
% dd__ = find(nf__==[ 200 500 1000 10000 20000 ]) ; 

% % ! sed -i '' -e "s/`head -9 Main.m | tail -1 `/nnx=110/" Main.m
% % system('sed -e "9s/.*/31419/" < Main.m > Main.m')
% % ! sed -i '' -e '9s/.*/replacement-line/' Main.m
% % 
% % eval(['!sed -i.bu ',char(39),'9s/.*/nn =   ',num2str(ii__),...
% %     '; /',char(39),' Main.m']);
% % 
% % ! sed -i.bu 's/oldword/newword/' file1.txt

Main ; 


% eval(['! cp out/RR_.log out_results/RR_',num2str(ii__),'_',num2str(dd__),'.log']);

eval(['! cp -r out out_results/out',num2str(ii__),'_',num2str(nf__),'']);


% pause
end

end

end