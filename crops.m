%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang, 2020
% The version of Matlab for this code is R2016b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = crops(in,N0,M)

out = zeros(N0,N0);

for ii=1:N0
    for jj=1:N0
        out(ii,jj) = in(ii+floor((M-N0)/2),jj+floor((M-N0)/2));%Cropped hologram
    end
end 
