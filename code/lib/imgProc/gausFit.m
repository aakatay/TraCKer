

%clear all; close all;
%I=double(imread('gaus.tif'));%assume gray scale, not RGB

function [cc_,RESNORM_,RESIDUAL_,EXITFLAG_,x] = gausFit(varargin)
    if nargin == 5
        I = cell2mat(varargin(1));
        fun = cell2mat(varargin(2));
        c0 = cell2mat(varargin(3));
        lb = cell2mat(varargin(4));
        ub = cell2mat(varargin(5));
    elseif nargin == 3
        I = cell2mat(varargin(1));
        fun = cell2mat(varargin(2));
        c0 = cell2mat(varargin(3));
        lb = [];
        ub = [];
    end

    I = double(I);
    %lb= [];ub=[];
    tolX = 0.2; % pixels
    tolX = 1e-6; % pixels
    TolFun = 1e-6; % pixels
    

    [m,n]=size(I);%assumes that I is a nxm matrix
    [X,Y]=meshgrid(1:n,1:m);%your x-y coordinates
    x(:,1)=X(:); % x= first column
    x(:,2)=Y(:); % y= second column
    f=I(:); % your data f(x,y) (in column vector)
    %--- now define the function in terms of x
    %--- where you use x(:,1)=X and x(:,2)=Y

    %--- now solve with lsqcurvefit
    %options = optimset('TolX',tolX,'Display','off','Algorithm','trust-region-reflective');  % 'Algorithm','levenberg-marquardt'
    options = optimset('TolX',tolX,'TolFun',TolFun,'Display','off','Algorithm','trust-region-reflective');
    %cc=lsqcurvefit(fun,c0,x,f);
    [cc, RESNORM,RESIDUAL,EXITFLAG]=lsqcurvefit(fun,c0,x,f,lb,ub,options);
    cc_ = cc;
    RESNORM_ = RESNORM;
    RESIDUAL_ = RESIDUAL;
    EXITFLAG_ = EXITFLAG;
end

