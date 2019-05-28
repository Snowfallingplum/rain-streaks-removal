%--------------Brief description-------------------------------------------
% This demo contains the implementation of the algorithm for video rain streaks removal
clear all;close all;clc;
path(path,genpath(pwd));
%%%--- Load Video ---%%%
%Assume that the number of rows and columns in the picture are equal
frames = 30;
load girl_rainy_45_heavy.mat
Rainy=Rainy(1:240,1:240,:,1:frames);
load girl_clean.mat
B_clean = B_clean(1:240,1:240,:,1:frames);
methodname{1}='   Rainy      ';
methodname{2}='DTV';
padsize=5;
%%
%implay(Rainy);implay(B_clean);
[O_Rainy,~]=rgb2gray_hsv(Rainy);%rgb2hsv
[O_clean,O_hsv]=rgb2gray_hsv(B_clean);
Rain=O_Rainy-O_clean;
PSNR0=PSNR3D(O_Rainy,O_clean);

SSIM_B10=ssim2(O_Rainy,O_clean);
SSIM_R10=ssim2(zeros(size(Rain)),Rain);
RSE0=norm(O_Rainy(:)-O_clean(:),'fro');

%% estimated direction
Dt = def3Dz;
direction=est_direction_patch(Dt(O_Rainy));
disp(direction);
%%%--- Parameters ---%%%
opts.tol=1e-2;
w_weight=[1 1 1];
opts.weight=w_weight/sum(w_weight);

opts.alpha1=1000;%DdR
opts.alpha2=100;%R
opts.alpha3=10;%DxB
opts.alpha4=100;%T
opts.alpha5=1;%L
opts.alpha6=opts.alpha3;%DyB
opts.tol= 1e-2;
opts.beta=50;
opts.maxit=70;
%---  ---%
O_Rainy = biger(O_Rainy,padsize);
%--- Rain streaks removal ---%
tic
[B_1,~]=rain_removal(O_Rainy,opts,direction);
time = toc;
%---  ---%
B_1 = smaller(B_1,padsize);
O_Rainy=smaller(O_Rainy,padsize);
R_1=O_Rainy-B_1;

%% reporting PSNR RSE SSIM of B and Rain
PSNR1= PSNR3D(B_1,O_clean); %psnr_all(kkk)=PSNR1;
SSIM_B11=    ssim2(B_1,O_clean); %ssim_all(kkk)=SSIM_B11;
SSIM_R11=    ssim2(R_1,Rain);
RSE1=    norm(B_1(:)-O_clean(:),'fro');
fprintf('\n');
fprintf('===========Time:  %5.3f=========================\n',time);
fprintf('        ||    %6.9s   ||%6.6s  ||  %6.7s  ||  %6.7s  ||  %6.6s        ||\n',' item ',' PSNR ', 'SSIM-B','SSIM-R',' RSE ');
fprintf('        || %6.9s  || %5.3f || %5.6f || %5.6f || %5.6f ||\n',...
    methodname{1},...
    PSNR0,...
    SSIM_B10,...
    SSIM_R10,...
    RSE0);
fprintf('        || %6.9s || %5.3f || %5.6f || %5.6f || %5.6f ||\n',...
    methodname{2},...
    PSNR1,...
    SSIM_B11,...
    SSIM_R11,...
    RSE1);
fprintf('===================================================\n');
%%%%
B_c=gray2color_hsv(O_hsv,B_1);
implay(B_c);
