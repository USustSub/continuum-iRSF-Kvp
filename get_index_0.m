 dd=[ij-1 ij ij+1 ij-1+nnx ij-1+nnx+1 ij-1+nnx+2  ij-1-nnx ij-1-nnx+1 ij-1-nnx+2];
            
            if ii == 1 
                index2 = sort(dd([1 2 4 5 7 8 ]));
%                 if norm(index-index2)>1e-10
%                     error('errr')
%                 end
%             Bp = B;
%             [~,B,en] = get_data ( xx , node , di , form ) ;
%             if norm(B-Bp)>1e-10
%                 error('eee')
%             end
            end
            
            if ii == 2 
                index2 = sort(dd([2 3  5 6  8 9]));
%                 if norm(index-index2)>1e-10
%                     error('errr')
%                 end
            end
            
            if ii == 3
                index2 = sort(dd([1 2 3 7 8 9  ]));
%                 if norm(index-index2)>1e-10
%                     error('errr')
%                 end
            end

            if ii == 4
                index2 = sort(dd([1 2 3 4 5 6  ]));
%                 if norm(index-index2)>1e-10
%                     error('errr')
%                 end
            end
            index = index2  ; 