%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang, 2020
% The version of Matlab for this code is R2016b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P] = Propagator_function(N, wavelength, area, z)

P = zeros(N,N);

for i = 1:N
    for j = 1:N
        alpha = wavelength*(i - N/2 -1)/area;
        beta = wavelength*(j - N/2 -1)/area;
        if ((alpha^2 + beta^2) <= 1)
        P(i,j) = exp(-2*pi*1i*z*sqrt(1 - alpha^2 - beta^2)/wavelength);
        end
    end
end
