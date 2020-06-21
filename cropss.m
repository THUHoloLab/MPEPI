%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang, 2020
% The version of Matlab for this code is R2016b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = cropss(in,N1,N2,M1,M2)
%N0 = 600;                  % initial number of pixels
%=2000;
out = zeros(M1-N1,M2-N2);
    
    for ii=1:(M1-N1)
        for jj=1:(M2-N2)
            out(ii,jj) = in(ii+N1,jj+N2);
        end
    end 
    
    
    
