%% Loads a data set with B_3 magnetic data and corresponding coordinates
%% Subseqeuently calculates the 3D-Hilbert transform both in frequency and space domain
%% to get B_1 and B_2 from B_3. The space domain approach is solved with the 
%% fastsum algorithm.

%Add libraries and clear workspace:
close all;
clear all;
addpath('..\nfft-3.5.1-mexw64-openmp\fastsum');
addpath('..\Util');
addpath('.\AntonSemechko-Bounding-Spheres-And-Circles-9555fec');

%%%%%%%%%%%%%%%
%% Paramters %%
%%%%%%%%%%%%%%%
load('sampledata_Kirchberg.mat');
rSurface = [RW,HW];

%Remove duplicates
[rSurface, ind] = unique(rSurface, 'rows');
deltaT = deltaT(ind);
deltaZ = deltaZ(ind);

N = length(deltaT);

year = 1965;
[a0, b0, c0, X0, Y0, Z0] = getTotalfield(year);

eps_B = 1/16;   % outer boundary
R = 0.25-eps_B/2; %actual radius of circle
[R0,C0,Xb] = ExactMinBoundCircle(rSurface);

%Move and Scale coordinates
rSurfaceScaled = (rSurface-C0)/R0*R;

%%%%%%%%%%%%%%%%%%%%%%
%% Frequency domain %%
%%%%%%%%%%%%%%%%%%%%%%
timeFreq = 0;
tf = tic;
xFreq = linspace(min(RW), max(RW), round(sqrt(N)));
yFreq = linspace(min(HW), max(HW), round(sqrt(N)));
[xmesh,ymesh] = meshgrid(xFreq,yFreq);
xyFreq = [xmesh(:),ymesh(:)];
dxFreq = min(uniquetol(diff(xFreq)));
dyFreq = min(uniquetol(diff(yFreq)));
nxFreq = length(xFreq);
nyFreq = length(yFreq);

F = scatteredInterpolant(rSurface(:,1),rSurface(:,2),deltaZ);
F.Method = 'linear';
B3Freq = F(xyFreq(:,1),xyFreq(:,2));   

[KX, KY, KXY] = calculateWavenumbers(nxFreq, nyFreq, dxFreq, dyFreq);
[HX, HY] = provideHilbertOperator(KX, KY);
[dY, dX] = performHilbertTransform(reshape(B3Freq,nxFreq,nyFreq), HX, HY);

F = scatteredInterpolant(xyFreq(:,1),xyFreq(:,2),dX(:));
F.Method = 'linear';
B1Freq = F(rSurface(:,1),rSurface(:,2));   
F = scatteredInterpolant(xyFreq(:,1),xyFreq(:,2),dY(:));
F.Method = 'linear';
B2Freq = F(rSurface(:,1),rSurface(:,2));   
BtFreq = B1Freq * a0 + B2Freq * b0 + deltaZ * c0;
timeFreq = toc(tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Space domain (fastsum) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% non-uniform source/target nodes on circle 
xFastsum = rSurfaceScaled;

% Calculate target nodes and Weighting of source nodes
[A, yFastsum] = getVoronoiWeightingAndRefPoints(xFastsum, 10, 'off', 0.25-eps_B/2);
M = length(yFastsum(:,1));


boolFast = 1; %set 1 for fastsum, set 0 for direct calculation
%Perform fastsum
B1Fast = hiFastsum(xFastsum, yFastsum, 1, eps_B, deltaZ, A, boolFast);
B2Fast = hiFastsum(xFastsum, yFastsum, 2, eps_B, deltaZ, A, boolFast);
BtFast = B1Fast * a0 + B2Fast * b0 + deltaZ * c0;
% Error
absErr =  abs(BtFast-deltaT);
relErr = absErr./abs(deltaT);
tFast = toc;
 

%%%%%%%%%%%%%
%% Figures %%
%%%%%%%%%%%%%
%define parameter to plot
pAna = deltaT;
pCalc = BtFast;
pError = abs(pAna-pCalc);

disp('Start plotting');
close all;
clf('reset');

mini = -300;
maxi = 1000;

% define plots
%parameters:
loadPlotParameters;

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
axis tight;
view(0, 90);
hold on;
scatter3(xFastsum(:,1),xFastsum(:,2),pAna,5,pAna,'filled');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
set(gca, 'CLim', [mini maxi]);
colormap('bluewhitered');
title('Bt (frequency)');

figure(2);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
axis tight;
view(0, 90);
hold on;
scatter3(xFastsum(:,1),xFastsum(:,2),pCalc,5,pCalc,'filled');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
set(gca, 'CLim', [mini maxi]);
colormap('bluewhitered');
title('Bt (fastsum)');


figure(3);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
scatter3(xFastsum(:,1),xFastsum(:,2),pError,5,pError,'filled');
axis tight;
view(0, 90);
hold on;
set(gca, 'CLim', [0 100]);
colormap('bluewhitered');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
title('Abs. Error (freq-fast)');
hold off;


figure(4)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
scatter3(xFastsum(:,1),xFastsum(:,2),deltaZ,5,deltaZ,'filled');
axis tight;
view(0, 90);
hold on;
h = colorbar();
set(gca, 'CLim', [mini maxi]);
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
colormap('bluewhitered');
title('B3 (data)');
hold off;

