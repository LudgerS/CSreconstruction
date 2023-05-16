function [reco, lambdaC, fC, N, memory] = CSreco_automatic(kspData, sigma, lambdaStart, wPhi, wTV, nMaxIllinois, TOL, epsFrac, minLambda)
% Yields a compressed sensing reconstruction with automatic optimization
% of the regulatization strength based on the noise level following the 
% Discrepancy Principle.
%
% The Discrepancy Principle is based on the assumption that the optimal
% reconstruction should deviate from the measured data approximately by the
% same amount as the data deviates from the true, but unknown signal due to
% noise. Multiple reconstructions are computed adjusting the regularization
% strength each time to acheive the desired deviation. After finding a
% valid starting interval (too much regularization, too little
% regularization), the Illinois Algorithm (a version of Regula Falsi) is
% used for fine tuning.
%
% Morozov VA. On the solution of functional equations by the method of 
% regularization. Paper presented at: Dokl Akad Nauk SSSR; 1966; Moscow, USSR
%
% Kilmer ME, O’Leary DP. Choosing regularization parameters in iterative 
% methods for ill-posed problems. SIAM J Matrix Anal Appl. 2001;22:1204–1221
%
% Dowell M, Jarratt P. A modified regula falsi method for computing
% the root of an equation. BIT Numer Math. 1971;11:168–174.
%
% Argumemnts:
%
%   kspData - k-space data with low frequencies in center and no leading
%   empty dimensions (suqeeze(kspData) == kspData).
%
%   sigma - noise level of the k-space data, standard deviation of the
%   Gaussian noise in each (!) the real and imaginary channel
%
%   lambdaStart - starting guess for the regulatization strength. 
%   Reasonable values depend on the noise level. Typical values would be 
%   between 10^-5 and 10^-1. A good starting value is crucial for fast
%   convergence.
%
%   wPhi and wTV govern the weighting of the image l1-norm and isotropic
%   Total Variation regulatization. wPhi = 1, wTV = 0 -> only image
%   l1-norm; wPhi = 0, wTV = 1 -> only TV; wPhi = 1, wTV = 1 -> equal
%   weighting
%   For a constant c, a reconstruction with [c*lambda, wPhi, wTV] would
%   yield identical results as one with [lambda, c*wPhi, c*wTV].
%
%   nMaxIllinois - maximum number of reconstructions performed for the 
%   regulatization strength optimization
%   
%   TOL - relative tolerance in the optimization
%
%   epsFrac -  desired deviation of the reconstruction from the measured 
%   data as a fraction of the noise level. Recommended value = .97
%
%   minLambda - for very high SNR data the allowed deviation becomes very
%   small. Limiting the minimal strength of the regularization ensures
%   successful removal of aliasing artifacts
%

%% normalization
Max = max(abs(kspData(:)));
kspData = kspData/Max;
sigma = sigma/Max;


%% change of Fourier transform convention, find undersampling mask
kspData = bfftn(cifftn(kspData));   % input is centered, output follows Matlab's convention
mask = abs(kspData) > 10^(-12);

fprintf('Acceleration %.2f\n', numel(mask)/sum(mask(:)))


%% determine data dimension
dim = size(kspData);

if length(dim) == 2
    dimSwitch = '2D';
elseif length(dim) == 3
    dimSwitch = '3D';
else
    error('unsupported dimension')
end


%% parameters
m = 0.5;        % controls Illinois update, 0.5 assumes symmetry
mu = 1;         % ADMM step size
admmIter = 50;   


%% Precalculation
epsNull = 2*sum(mask(:))*sigma^2;                                               % expected deviation of the data from the true signal
devFunc = @(reco) sumN(abs(kspData - mask.*bfftn(reco)).^2) - epsFrac*epsNull;  % (deviation between reco and data) - (desired deviation)

switch dimSwitch
    case '2D'
        laplaceColumn = laplaceColumn2DnonIso(dim(1), dim(2));

    case '3D'
        laplaceColumn = laplaceColumn3DnonIso(dim(1), dim(2), dim(3));

    otherwise
        error('unsupported dimension')
end

eVs = sqrt(prod(dim))*bifftn(reshape(laplaceColumn, dim));
inverseO = 2*mask + wPhi^2*1 + wTV^2*eVs; 
inverseO = inverseO.^-1;


%% Reco function
switch dimSwitch
    case '2D'
        recoFunc = @(lambda) ADMM2D(kspData,inverseO,admmIter,mu,lambda,wPhi,wTV);

    case '3D'
        recoFunc = @(lambda) ADMM3D(kspData,inverseO,admmIter,mu,lambda,wPhi,wTV);

    otherwise
        error('unsupported dimension')
end


%% Find starting interval
N = 1;
fprintf('%d\n',N)

lambdaA = max(lambdaStart, minLambda);
reco = recoFunc(lambdaA);
fA = devFunc(reco);
fprintf('+      %.2e || %.4f\n',lambdaA,fA/epsNull)
Switch = 0;

memory = zeros(nMaxIllinois,2);
memory(N,:) = [lambdaA,fA];

if abs(fA) < (TOL*epsNull)
    fC = fA;
    lambdaC = lambdaA;
    prepareReturn('achieved')
    return
end

if fA < 0   % search upper border
    
    lambdaB = lambdaA;
    
    while (Switch == 0) && (N < nMaxIllinois)

        N = N+1;
        fprintf('%d\n',N)

        % new reco
        lambdaB = 2*lambdaB;
        reco = recoFunc(lambdaB);
        fB = devFunc(reco);

        % update on screen and saved information
        fprintf('+      %.2e || %.4f\n', lambdaB, fB/epsNull)
        memory(N,:) = [lambdaB, fB];

        % test deviation from target
        if abs(fB) < (TOL*epsNull)
            fC = fB;
            lambdaC = lambdaB;
            prepareReturn('achieved')
            return
        end

        % test if upper border was found
        if fB > 0
            Switch = 1;
        else
            lambdaA = lambdaB;
            fA = fB;
        end
        
    end
    
elseif fA > 0   % search lower border

    lambdaB = lambdaA;
    fB = fA;

    while (Switch == 0) && (N < nMaxIllinois) && (lambdaA > minLambda)

        N = N+1;
        fprintf('%d\n',N)

        % new reco
        lambdaA = max(0.5*lambdaA, minLambda);
        reco = recoFunc(lambdaA);
        fA = devFunc(reco);

        % update on screen and saved information
        fprintf('+      %.2e || %.4f\n',lambdaA,fA/epsNull)
        memory(N,:) = [lambdaA,fA];

        % test deviation from target
        if abs(fA) < (TOL*epsNull)
            fC = fA;
            lambdaC = lambdaA;
            prepareReturn('achieved')
            return
        end

        % test if lower border was found
        if fA < 0
            Switch = 1;
        else
            lambdaB = lambdaA;
            fB = fA;
        end
        
    end
    
    % return minLambda reco
    if lambdaA == minLambda
        fC = fA;
        lambdaC = lambdaA;
        prepareReturn('minLambda')            
        return
    end
    
end

if Switch == 0

    fprintf('\n Valid starting interval not found in %d trials \n',N)
    reco = NaN;
    lambdaC = NaN;
    fC = NaN;
    N = NaN;
    memory = NaN;

else
    
    
%% Illinois algorithm
N = N + 1;
fprintf('%d\n',N) 

% new reco
lambdaC = (fB*lambdaA - fA*lambdaB)/(fB - fA);
reco = recoFunc(lambdaC);
fC = devFunc(reco);

% update on screen and saved information
fprintf('+      %.2e || %.4f\n',lambdaC,fC/epsNull)
memory(N,:) = [lambdaC,fC];

while N < nMaxIllinois
    
    if abs(fC) < (TOL*epsNull)
        
        prepareReturn('achieved')   
        return

    elseif sign(fA) == sign(fC)     % last C becomes the new lower border

        lambdaA = lambdaC;
        fA = fC;
        lambdaC = (m*fB*lambdaA - fA*lambdaB)/(m*fB-fA);

    elseif sign(fB) == sign(fC)     % last C becomes the new upper border

        lambdaB = lambdaC;
        fB = fC;
        lambdaC = (fB*lambdaA - m*fA*lambdaB)/(fB-m*fA);

    end

    N = N + 1;
    fprintf('%d\n',N) 
    
    reco = recoFunc(lambdaC);
    fC = devFunc(reco);
    fprintf('+      %.2e || %.4f\n', lambdaC, fC/epsNull)
      
    memory(N,:) = [lambdaC,fC];

end

end

% return last reconstruction
prepareReturn('notAchieved')


function [] = prepareReturn(switchValue)
% clean up, rescale and print appropriate on-screen message
    switch switchValue
        case 'achieved'
            fprintf('tolerance achieved\n')
        case 'minLambda'
            fprintf('min lambda reached\n')
        case 'notAchieved'
            fprintf('tolerance not achieved\n')
            fprintf('check memory and consider increasing nMaxIllinois\n')
        otherwise
            error('wrong switch statement')
    end
            
    memory = memory(1:N,:);
    reco = reco*Max;
    fprintf('\n')
end





end

