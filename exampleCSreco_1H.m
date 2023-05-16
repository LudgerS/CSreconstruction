%% exampleCSreco_1H.m
%
% part of https://github.com/LudgerS/CSreconstruction  
%
% Example script computing compressed sensing reconstructions of the 1H
% MRI example data with automatic optimization and manual setting of the
% regularization strength.
%
% To run the script, ensure that exampleData.mat and the 'subFunctions' 
% folder are on the search path. This should automatically be the case if 
% the repository is simply downloaded and the scripts location is made the 
% curent folder.
% 
% For details on the algorithm, please consider the readme and 
%   Starke, Ludger, et al. 
%   "Performance of compressed sensing for fluorine?19 magnetic resonance 
%   imaging at low signal?to?noise ratio conditions." 
%   Magnetic resonance in medicine 84.2 (2020): 592-608.
%
% Written by Ludger Starke
%
% License: GNU GPLv3 


clear, close all
addpath(genpath([pwd, filesep, 'subFunctions']))


%% load example data and select nucleus
exampleData = load('exampleData.mat');

ksp = exampleData.ksp1H;       % ksp-data
sigma = exampleData.sigma1H;   % noise level

dim = size(ksp);


%% create and apply mask for retrospective undersampling
usFactor = 1/3;      % inverse of acceleration

centerFraction = 0.1;   % mask params. defined as in paper, can be kept standard
degree = 1.5;

usMask_phasePlane = polynomial2DsamplingPattern(dim(2:3), usFactor, centerFraction, degree);
usMask = permute(repmat(usMask_phasePlane, [1, 1, dim(1)]), [3, 1, 2]);

ksp_us = ksp;
ksp_us(~usMask) = 0;


%% conventional Fourier reconstruction
fullySampled = cifftn(ksp);
zeroFilled = cifftn(ksp_us);


%% CS reconstruction with automatic lambda
lambdaStart = 3.5*10^-5;    % starting value of regularizaton strength optimization
minLambda = 10^(-7);        % minimal allowed lambda value

% weighting of image l1-norm (wPhi) and Total Variation (wTV) regularizaton
wPhi = 0;   %  The 1H image is not sparse without a transformation
wTV = 1;

nMaxIllinois = 20;      % maximum number of steps used to optimize lambda 
epsilonFrac = 0.97;     % desired deviation of the reconstruction from the 
                        % measured data (epsilon) as a fraction of the noise level
TOL = 0.01;             % relative tolerance regarding epsilon

[CSreco_auto, lambdaC, fC, N, memory] = CSreco_automatic(ksp_us, sigma, lambdaStart, wPhi, wTV, nMaxIllinois, TOL, epsilonFrac, minLambda);


%% CS reconstruction with manual lambda
lambda = 6*10^-5;
CSreco_man = CSreco_manual(ksp_us, lambda, wPhi, wTV);


%% show result
concatenatedRecos = cat(2, absNorm(fullySampled), absNorm(zeroFilled), absNorm(CSreco_auto), absNorm(CSreco_man));
titleStrings = ["fully-sampled", "zero-filled", "CS auto", "CS manual"];
display3dImage(concatenatedRecos, gray(512), 3, titleStrings)

