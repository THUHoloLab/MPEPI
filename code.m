%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang, 2020
% The version of Matlab for this code is R2016b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all;
clear all;

N1 = 400;                  % initial number of pixels
N = 1400;
N2 = 1400;

object=objectPro0(N);
phase=phasePro0(N);
object=object.*(exp(1i.*phase));
support=supportPro0(N2);

wavelength = 500*10^(-9); 
pixel=1*10^(-6); %Pixel size
area = N*pixel; %Cover area of initial whole hologram
z0 = 0.0022;%Plane5
z1 = 0.0021;%Plane4
z2 = 0.0020;%Plane3
z3 = 0.0019;%Plane2
z4 = 0.0018;%Plane1

prop0 = Propagator_function(N, wavelength, area, z0);
prop1 = Propagator_function(N, wavelength, area, z1);
prop2 = Propagator_function(N, wavelength, area, z2);
prop3 = Propagator_function(N, wavelength, area, z3);
prop4 = Propagator_function(N, wavelength, area, z4);

area1 = N1*pixel; %Cover area of small size hologram
prop00 = Propagator_function(N1, wavelength, area1, z0);
prop11 = Propagator_function(N1, wavelength, area1, z1);
prop22 = Propagator_function(N1, wavelength, area1, z2);
prop33 = Propagator_function(N1, wavelength, area1, z3);
prop44 = Propagator_function(N1, wavelength, area1, z4);

area2 = N2*pixel; %Cover area of final size hologram
prop000 = Propagator_function(N2, wavelength, area2, z0);
prop111 = Propagator_function(N2, wavelength, area2, z1);
prop222 = Propagator_function(N2, wavelength, area2, z2);
prop333 = Propagator_function(N2, wavelength, area2, z3);
prop444 = Propagator_function(N2, wavelength, area2, z4);
U0 = IFT(FT(object).*prop0);
U1 = IFT(FT(object).*prop1);
U2 = IFT(FT(object).*prop2);
U3 = IFT(FT(object).*prop3);
U4 = IFT(FT(object).*prop4);
hologram0 = abs(U0).^2;%hologram
hologram1 = abs(U1).^2;
hologram2 = abs(U2).^2;
hologram3 = abs(U3).^2;
hologram4 = abs(U4).^2;
figure(3),imshow(abs(hologram0),[])


%%

hologram_cropped0=crops(hologram0,N1,N);% Capture small size hologram
hologram_cropped1=crops(hologram1,N1,N);
hologram_cropped2=crops(hologram2,N1,N);
hologram_cropped3=crops(hologram3,N1,N);
hologram_cropped4=crops(hologram4,N1,N);
recons0 = IFT(FT(sqrt(hologram_cropped0)).*conj(prop00));%Reconstruction of small size hologram
recons1 = IFT(FT(sqrt(hologram_cropped1)).*conj(prop11));
recons2 = IFT(FT(sqrt(hologram_cropped2)).*conj(prop22));
recons3 = IFT(FT(sqrt(hologram_cropped3)).*conj(prop33));
recons4 = IFT(FT(sqrt(hologram_cropped4)).*conj(prop44));
rec0=crops(recons0,50,N1);
rec1=crops(recons1,50,N1);
rec2=crops(recons2,50,N1);
rec3=crops(recons3,50,N1);
rec4=crops(recons4,50,N1);
figure(3),subplot('position',[0.01 0.01 0.98 0.98]),imshow(abs(hologram_cropped4),[])

measured0=sqrt(hologram_cropped0);% small size hologram
measured1=sqrt(hologram_cropped1);
measured2=sqrt(hologram_cropped2);
measured3=sqrt(hologram_cropped3);
measured4=sqrt(hologram_cropped4);
amplitude0 = ones(N2,N2);phase0 = zeros(N2,N2);
amplitude1 = ones(N2,N2);phase1 = zeros(N2,N2);
amplitude2 = ones(N2,N2);phase2 = zeros(N2,N2);
amplitude3 = ones(N2,N2);phase3 = zeros(N2,N2);
amplitude4 = ones(N2,N2);phase4 = zeros(N2,N2);
field_final0 = amplitude0.*exp(1i.*phase0);
field_final1 = amplitude1.*exp(1i.*phase1);
field_final2 = amplitude2.*exp(1i.*phase2);
field_final3 = amplitude3.*exp(1i.*phase3);
field_final4 = amplitude4.*exp(1i.*phase4);

M1 = (N2 - N1)/2;
M2 = (N2 + N1)/2;

%% Iteration
tic
Iterations=400;% Iterative number
for kk = 1:Iterations

fprintf('Iteration: %d\n', kk)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plane5
for ii = M1+1:M2
for jj = M1+1:M2
  amplitude0(ii,jj) = measured0(ii-M1,jj-M1);

end
end
field_final0 = amplitude0.*exp(1i.*phase0);
t0 = IFT((FT(field_final0)).*conj(prop000));

object0 = 1 - abs(t0); % Object amplitude
ph0=angle(t0);% Object phase
for ii=1:N2
    for jj=1:N2
        if (object0(ii,jj)<0)
            object0(ii,jj)=0;
            ph0(ii,jj)=0;
        end
    end
end
object0 = object0.*support;
t0 = (1 - object0).*exp(1i*ph0);

field_final_updated1 = IFT((FT(t0)).*(prop111));% Propagate to the next hologram plane
amplitude1 = abs(field_final_updated1);% Undated amplitude in next hologram plane
phase1 = angle(field_final_updated1);% Undated phase in next hologram plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plane4
for ii = M1+1:M2
for jj = M1+1:M2
  amplitude1(ii,jj) = measured1(ii-M1,jj-M1);

end
end
field_final1 = amplitude1.*exp(1i.*phase1);
t1 = IFT((FT(field_final1)).*conj(prop111));

object1 = 1 - abs(t1); 
ph1=angle(t1);
for ii=1:N2
    for jj=1:N2
        if (object1(ii,jj)<0)
            object1(ii,jj)=0;
            ph1(ii,jj)=0;
        end
    end
end
object1 = object1.*support;
t1 = (1 - object1).*exp(1i*ph1);

field_final_updated2 = IFT((FT(t1)).*(prop222));% Propagate to the next hologram plane
amplitude2 = abs(field_final_updated2);% Undated amplitude in next hologram plane
phase2 = angle(field_final_updated2);% Undated phase in next hologram plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plane3
for ii = M1+1:M2
for jj = M1+1:M2
  amplitude2(ii,jj) = measured2(ii-M1,jj-M1);

end
end
field_final2 = amplitude2.*exp(1i.*phase2);
t2 = IFT((FT(field_final2)).*conj(prop222));

object2 = 1 - abs(t2); 
ph2=angle(t2);
for ii=1:N2
    for jj=1:N2
        if (object2(ii,jj)<0)
            object2(ii,jj)=0;
            ph2(ii,jj)=0;
        end
    end
end
object2 = object2.*support;
t2 = (1 - object2).*exp(1i*ph2);

field_final_updated3 = IFT((FT(t2).*(prop333)));% Propagate to the next hologram plane
amplitude3 = abs(field_final_updated3);% Undated amplitude in next hologram plane
phase3 = angle(field_final_updated3);% Undated phase in next hologram plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plane2
for ii = M1+1:M2
for jj = M1+1:M2
  amplitude3(ii,jj) = measured3(ii-M1,jj-M1);

end
end
field_final3 = amplitude3.*exp(1i.*phase3);
t3 = IFT((FT(field_final3)).*conj(prop333));

object3 = 1 - abs(t3); 
ph3=angle(t3);
for ii=1:N2
    for jj=1:N2
        if (object3(ii,jj)<0)
            object3(ii,jj)=0;
            ph3(ii,jj)=0;
        end
    end
end
object3 = object3.*support;
t3 = (1 - object3).*exp(1i*ph3);

field_final_updated4 = IFT((FT(t3).*(prop444)));% Propagate to the next hologram plane
amplitude4 = abs(field_final_updated4);% Undated amplitude in next hologram plane
phase4 = angle(field_final_updated4);% Undated phase in next hologram plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plane1
for ii = M1+1:M2
for jj = M1+1:M2
  amplitude4(ii,jj) = measured4(ii-M1,jj-M1);

end
end
field_final4 = amplitude4.*exp(1i.*phase4);
t4 = IFT((FT(field_final4)).*conj(prop444));

object4 = 1 - abs(t4); 
ph4=angle(t4);
for ii=1:N2
    for jj=1:N2
        if (object4(ii,jj)<0)
            object4(ii,jj)=0;
            ph4(ii,jj)=0;
        end
    end
end
object4 = object4.*support;
t4 = (1 - object4).*exp(1i*ph4);

field_final_updated0 = IFT((FT(t4).*(prop000)));% Propagate to the field_final0 plane
amplitude0 = abs(field_final_updated0);% Undated amplitude in the field_final0 plane
phase0 = angle(field_final_updated0);% Undated phase in the field_final0 plane

if kk==200;
    tt0=t0;
end

end
toc
t00=crops(t0,50,N2);
figure(3),subplot('position',[0.01 0.01 0.8 0.98]),imshow(angle(t00),[])
colorbar('location','EastOutside','FontSize',14)


