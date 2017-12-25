close all;
clear all

%% angular distribution
% parameters
n = 3; % # of diffraction orders
sq = 4; % squeeze factor
wl = 0.5*1e-6;
szPinhole = 50*1e-6; % pinhole diameter
frstzero = 1.22*wl/(szPinhole); % [radians]
offset = 0; % set hard edge coordinate
nTheta = 1e4; % sampling

k = 2*pi/wl; 
r = szPinhole/2; % pinhole radius

% profile
t = linspace(0,frstzero*n,nTheta); % theta [radians]
x_ = k*r*sin(t); % bessel x
% Gaussian
I1Dhalf = (2*besselj(1,x_)).^2./x_.^2;
% Squeezed Gaussian
I1DhalfSq = I1Dhalf(1:sq:nTheta);
I1DhalfSq(nTheta/sq/n-offset:end ) = 0;   
%combined halves
I1D = fliplr([fliplr(I1Dhalf) I1DhalfSq]);
%display
hFig = figure(1);
xx = -nTheta+1:nTheta/sq;
hp = plot(xx,I1D,'k','LineWidth',2,'LineStyle','-');

%% spacial distribution
% profile
t = linspace(0,frstzero*n,nTheta); % theta [radians]
x_ = k*r*sin(t); % bessel x
% Gaussian
I1Dhalf = (2*besselj(1,x_)).^2./x_.^2;
% Squeezed Gaussian
I1DhalfSq = I1Dhalf(1:sq:nTheta);
%combined halves
I1D = fliplr([fliplr(I1Dhalf) I1Dhalf]);
%display
hFig = figure(2);
xx = -nTheta+1:nTheta;
hp = plot(xx,I1D,'k','LineWidth',2,'LineStyle','-');

