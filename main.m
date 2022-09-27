%% BEST RUN WITH MATLAB R2018b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was solely written by Jianjun Jiao in Lanzhou Jiaotong University.
% Correspondence address:Anning West Road of Anning District, Lanzhou, China.
% Contact:jiaojianjun@lzjtu.edu.cn
% Last change: January, 2022.
% Basically, you can run this code SEVERAL times to acquire the most desired result.
% It is welcomed to change the following parameters as you like to see what gonna happen.
%% Inputs:
% cluster_num - Number of clusters
% max_iter -  Maximum number of iterations
% density  -  The noise density
% Rg - The radius of Gaussian filtering window
% Ra - The radius of average filtering window
% ==============Parameters for region-level information================
% l - Side length of block
% S - Side length of region
% h - Attenuation of exponential function in Eqs. (6)-(7)
% sigma - Gaussian standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
clc;
clear all;
close all;
%%  Input image
f_uint8 = imread('3063.jpg'); 
[m, n, p] = size(f_uint8);
I= double(f_uint8);
%%  Parameters setting
%Parameters for clustering
max_iter=200;
cluster_num =3;
Rg=9;
Ra=17;
% Parameters for region-level information
l = 9;
S = 15;
h = 7; 
sigma=4;
%% Construct mixed noise
density=0.20;
I= I / 255;
I = imnoise(I, 'gaussian', 0, density);
I = imnoise(I, 'salt & pepper', density);
I = imnoise(I, 'speckle', density);
figure, imshow(I),title('noise image');
I = I * 255;
%%  Acquire region_level information
I_R = region_level_information(I, l, S, h, sigma);
% figure, imshow(uint8(I_R)),title('I_R');
% Nomalization
I_R = I_R / 255;
I = I / 255;
%%  Gaussian and mean smoothing image
mask = fspecial('gaussian', [Rg Rg],sigma);
mask1 = fspecial('average',Ra);
I_G = imfilter(I_R, mask, 'symmetric');  
I_A = imfilter(I_R, mask1, 'symmetric'); 
% figure, imshow(I_G),title('I_G');
% figure, imshow(I_A),title('I_A');
%% Calculate  weights  
[row, col, depth] = size(I);
n = row * col;
phi = 255 * abs(I_G - I_R); % Eq. (9) 
phi_padded = padarray(phi, [1 1], 'symmetric');
C_j = zeros(row, col,depth);
k=1;%
r=3;%
f_padded=phi_padded;
%%%%%%%%%%%%%%% Eq.(10)
for i=-k:k
    for j=-k:k 
        
        C_j =C_j +abs(f_padded(i + k + 1 : end + i - k, j + k + 1 : end + j - k,:)-imfilter(f_padded(i + k + 1 : end + i - k, j + k + 1 : end + j - k,:),  mask, 'replicate'))/r*r; 
    
    end
end
alpha=C_j;
beta=C_j/r*r;
Ph=beta+alpha; %Eq. (13)
Ph_min = min(min(min(Ph)));
Ph_max = max(max(max(Ph)));
q = (Ph - Ph_min) ./ (Ph_max  -Ph_min);%Eq. (14)
gamma=1-q; %Eq. (12)

Xi=(1./(alpha+alpha+beta)).*(alpha.*I_A+beta.*I_G+gamma.*I_R);% Eq. (8)
Xi_min = min(min(min(Xi)));
Xi_max = max(max(max(Xi)));
Xi = (Xi - Xi_min) ./ (Xi_max  -Xi_min);
figure, imshow(Xi,[]),title('The weighted image');
%% Calculate Eq. (21)
I=I*255;
Xi=Xi*255;
I_R=I_R*255;
difference =50*abs(log(I-Xi+ eps)) ;
%% The way to set adaptive coefficient
Lambda=difference;%Eq. (18)
Omega=1./difference;%Eq. (17)
%% RLFCM
[U,C]= RLFCM(I,Xi,cluster_num,max_iter,Omega,Lambda);

