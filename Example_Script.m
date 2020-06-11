%-----------------------------------------------------------
% Continuous Fourier Ring Correlation
%
% Exemple of use of the function FRCcont with and without the use of the
% lookup table. 
%
% Copyright (2018) Thanh-An Pham (thanh-an.pham@epfl.ch)
%                  Emmanuel Soubies (esoubies@gmail.com)
%-----------------------------------------------------------
% clear, close all;

%% Parameters
fov=6400;         % Field of view to generate points (nm)     
Nmol=1000;        % Number of molecules
stdMol=8;         % Standard deviation between the position of the molecules of the two sets
Nsamples = 200;   % Number of points in the FRC curve
fmax = 0.1;       % Maximal frequency. The FRC will be computed on [0,fmax]

%% Load data
% gt=rand(Nmol,3)*fov;
% loc=gt+randn(Nmol,3)*stdMol;

gt=rand(Nmol,3)*fov;

nDetectsPerPoint = 10;

nPtsActual = repmat(nDetectsPerPoint, [Nmol, 1]);
nPtsActual = poissrnd(nPtsActual);
pts = zeros(sum(nPtsActual(:)), 4); %[x, y, z, nPhotons, locPrecision]
ptsNow = 1;
for k = 1:size(nPtsActual, 1)

    precHere = stdMol.^2;
    
    pts(ptsNow:(ptsNow + nPtsActual(k)-1),:) = ...
        [gt(k,1) + precHere.*randn(nPtsActual(k), 1), ...
         gt(k,2) + precHere.*randn(nPtsActual(k), 1), ...
         gt(k,3) + precHere.*randn(nPtsActual(k), 1), ...
         repmat(precHere, [nPtsActual(k), 1])] ;
     
     ptsNow = ptsNow + nPtsActual(k);
    
end

%% Split
pts = pts(randperm(size(pts, 1)), :);
s1 = pts(1:2:end, :);
s2 = pts(2:2:end, :);

%% FRC computation
tic;[FRCcLook,freq] = FRCcont(s1,s2,fov,Nsamples,fmax,1);
disp(['Computation with Lookup table : ',num2str(toc),' s']);
tic;[FRCcNoLook,freq] = FRCcont(s1,s2,fov,Nsamples,fmax,0);
disp(['Computation without Lookup table : ',num2str(toc),' s']);

%% Plots
figure;hold all;
plot(freq,FRCcLook,'linewidth',1.5);
plot(freq,FRCcNoLook,'linewidth',1.5,'linestyle','--');grid
xlabel('f'); ylabel('FRC'); set(gca,'FontSize',14);
legend('Lookup Table','No Lookup Table');axis([0 0.1 -0.05 1]);

