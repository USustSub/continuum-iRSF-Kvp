 dd=[ij-1 ij ij+1 ij-1+nnx ij-1+nnx+1 ij-1+nnx+2  ij-1-nnx ij-1-nnx+1 ij-1-nnx+2];
            
            if ii == 1 
                index2 = sort(dd([1 2 4 5 7 8 ]));
                en__l = zeros(1,2*size(index2,2));  
                for m = 1 : size(index2,2)      
                    en__l(2*m-1) = 2*index2(m)-1;
                    en__l(2*m  ) = 2*index2(m)  ;
                end
                if isempty(B__l)
%                     [~,B__l,~] = get_data ( xx , node , di , form ) ;
                    [~,B__l,~] = get_data_2 ( xx , node , index2 , form , di) ;
                end
                B = B__l ; 
                en = en__l ; 
            end

            if ii == 2
                index2 = sort(dd([2 3  5 6  8 9]));
                en__r = zeros(1,2*size(index2,2));  
                for m = 1 : size(index2,2)      
                    en__r(2*m-1) = 2*index2(m)-1;
                    en__r(2*m  ) = 2*index2(m)  ;
                end
                if isempty(B__r)
%                     [~,B__r,~] = get_data ( xx , node , di , form ) ;
                    [~,B__r,~] = get_data_2 ( xx , node , index2 , form , di) ;                    
                end
                B = B__r ; 
                en = en__r ; 
            end
            
            if ii == 3
                index2 = sort(dd([1 2 3 7 8 9  ]));
                en__b = zeros(1,2*size(index2,2));  
                for m = 1 : size(index2,2)      
                    en__b(2*m-1) = 2*index2(m)-1;
                    en__b(2*m  ) = 2*index2(m)  ;
                end
                if isempty(B__b)
%                     [~,B__b,~] = get_data ( xx , node , di , form ) ;
                    [~,B__b,~] = get_data_2 ( xx , node , index2 , form , di) ;                    
                end
                B = B__b ; 
                en = en__b ; 
            end

            if ii == 4
                index2 = sort(dd([1 2 3 4 5 6  ]));
                en__t = zeros(1,2*size(index2,2));  
                for m = 1 : size(index2,2)      
                    en__t(2*m-1) = 2*index2(m)-1;
                    en__t(2*m  ) = 2*index2(m)  ;
                end

                if isempty(B__t)
%                     [~,B__t,~] = get_data ( xx , node , di , form ) ;
                    [~,B__t,~] = get_data_2 ( xx , node , index2 , form , di) ;                                        
                end
                B = B__t ; 
                en = en__t ; 
            end
            index = index2 ; 