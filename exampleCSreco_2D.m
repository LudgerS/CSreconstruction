%% exampleCSreco_2D.m
%
% part of https://github.com/LudgerS/CSreconstruction  
%
% Example script a computing compressed sensing MRI reconstructions of 
% simulated 2D example data with both automatic optimization and manual 
% setting of the regularization strength.
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


%% create digital phantom
exampleData = load('exampleData.mat');

simulatedImage = exampleData.digitalPhantom;
dim = size(simulatedImage);

pSNR = 20;                                   % desired peak SNR in fully-sampled reco
sigma = max(simulatedImage(:))/pSNR;
simulatedImage = simulatedImage + sigma*(randn(dim) + 1i*randn(dim));

ksp = cfftn(simulatedImage);                 % ksp-data


%% create and apply mask for retrospective undersampling
usFactor = 1/3;         % inverse of acceleration

centerFraction = 0.1;   % mask params. defined as in paper, can be kept standard
degree = 1.5;

usMask = polynomial1DsamplingPattern(dim, usFactor, centerFraction, degree);

ksp_us = ksp;
ksp_us(~usMask) = 0;


%% conventional Fourier reconstruction
fullySampled = cifftn(ksp);
zeroFilled = cifftn(ksp_us);


%% CS reconstruction with automatic lambda
lambdaStart = 0.015;    % starting value of regularizaton strength optimization
minLambda = 10^(-7);    % minimal allowed lambda value

% weighting of image l1-norm (wPhi) and Total Variation (wTV) regularizaton
wPhi = 1;
wTV = 1;

nMaxIllinois = 20;      % maximum number of steps used to optimize lambda 
epsilonFrac = 0.97;     % desired deviation of the reconstruction from the 
                        % measured data (epsilon) as a fraction of the noise level
TOL = 0.01;             % relative tolerance regarding epsilon

[CSreco_auto, lambdaC, fC, N, memory] = CSreco_automatic(ksp_us, sigma, lambdaStart, wPhi, wTV, nMaxIllinois, TOL, epsilonFrac, minLambda);


%% CS reconstruction with manual lambda
lambda = 0.017;
CSreco_man = CSreco_manual(ksp_us, lambda, wPhi, wTV);


%% show result
concatenatedRecos = cat(2, absNorm(fullySampled), absNorm(zeroFilled), absNorm(CSreco_auto), absNorm(CSreco_man));
titleStrings = ["fully-sampled", "zero-filled", "CS auto", "CS manual"];
display3dImage(concatenatedRecos, gray(512), 3, titleStrings)


