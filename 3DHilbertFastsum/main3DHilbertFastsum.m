%% Creates a model of a magnetized sphere and its magnetic field on the surface
%% Subseqeuently calculates the 3D-Hilbert transform both in frequency and space domain
%% to get B_1 and B_2 from B_3. The space domain approach is solved with the 
%% fastsum algorithm.

%Add libraries and clear workspace:
close all;
clear all;

addpath('..\nfft-3.5.1-mexw64-openmp\fastsum');
addpath('..\Util');


%%%%%%%%%%%
%% Model %%
%%%%%%%%%%%

%Parameters:
%Definde circle for surface points
N = 7200; %number of points
eps_B = 0.12;   % outer boundary
R = 0.25-eps_B/2; %actual radius of circle

% Non uniform surface/measurement points 
rSurface = randsphere(N,2,R,10);;
N = length(rSurface);



%% Subsurface Modell:
% The magnetic total field T0, given in units of nT at the time of the measurement, 
% is defined by its three cartesian components X0, Y0, and Z0.
% The components of the unit vector in the direction of T0 are indicated by a0, b0, and c0.
% We will consider a magnetic main field obtained from the IGRF valid for
% 1965.

% Now, let's define a magnetized sphere beneath the observation points.
% We consider the sphere to be located at x = ? [m], y = ß [m], and a depth 
% of t [m]. Its radius is a [m]. The magnetic susceptibility is kappa SI units.
t = 0.05;
a = 0.02;
rp = [0 0 t];
vol = 4 / 3 * pi * a^3;
kappa = 0.126;

year = 1965;
[a0, b0, c0, X0, Y0, Z0] = getTotalfield(year);
T0 = [X0 Y0 Z0];

% The total magnetization of the sphere is expressed by magn:
mu0 = 4e-7 * pi;
magn = T0 / mu0 * kappa * vol;

% Now, let's create data and store the cartesian components of the anomaly 
% in the array B.
B = zeros(3, N);
for i = 1:N
    B(:, i) = Bfield(magn, [rSurface(i,:) 0], rp);
end

% Finally, the total field anomaly can be obtained by projection of the 
% magnetic anomaly onto the total field.
Bt = B(1, :, :) * a0 + B(2, :, :) * b0 + B(3, :, :) * c0;


%%%%%%%%%%%%%%%%%%%%%%
%% Frequency domain %%
%%%%%%%%%%%%%%%%%%%%%%
timefft = 0;
tf = tic;
xFreq = linspace(-R, R, round(sqrt(N)));
yFreq = linspace(-R, R, round(sqrt(N)));
[xmesh,ymesh] = meshgrid(xFreq,yFreq);
xyFreq = [xmesh(:),ymesh(:)];
dxFreq = uniquetol(diff(xFreq));
dyFreq = uniquetol(diff(yFreq));
nxFreq = length(xFreq);
nyFreq = length(yFreq);

F = scatteredInterpolant(rSurface(:,1),rSurface(:,2),B(3, :, :)');
F.Method = 'linear';
B3Freq = F(xyFreq(:,1),xyFreq(:,2));   

[KX, KY, KXY] = calculateWavenumbers(nxFreq, nyFreq, dxFreq, dyFreq);
[HX, HY] = provideHilbertOperator(KX, KY);
[dY, dX] = performHilbertTransform(reshape(B3Freq,nxFreq,nyFreq), HX, HY);

F = scatteredInterpolant(xyFreq(:,1),xyFreq(:,2),dX(:));
F.Method = 'linear';
BxFreq = F(rSurface(:,1),rSurface(:,2));   
F = scatteredInterpolant(xyFreq(:,1),xyFreq(:,2),dY(:));
F.Method = 'linear';
ByFreq = F(rSurface(:,1),rSurface(:,2));   
BtFreq = BxFreq' * a0 + ByFreq' * b0 + B(3, :, :) * c0;
timefft = toc(tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Space domain (fastsum) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%non-uniform source/target nodes on circle 
xFastsum = rSurface;

% Calculate Weighting of source nodes
[A, y_fastsum] = getVoronoiWeightingAndRefPoints(xFastsum, 10, 'on', 0.25-eps_B/2);
M = length(y_fastsum(:,1));


boolFast = 1;
tstart = tic;
%Perform fastsum
BxFast = hiFastsum(xFastsum, y_fastsum, 1, eps_B, B(3, :, :)', A, boolFast);
ByFast = hiFastsum(xFastsum, y_fastsum, 2, eps_B, B(3, :, :)', A, boolFast);
BtFast = BxFast * a0 + ByFast * b0 + B(3, :, :)' * c0;

tend = toc;
% Error
absErr =  abs(BtFast-Bt');
relErr = absErr./abs(Bt');

%%%%%%%%%%%%%
%% Figures %%
%%%%%%%%%%%%%
%define parameter to plot
pCalc = BtFast;
pAna = Bt';
pError = abs(pAna-pCalc);

disp('Start plotting');
close all;
clf('reset');

mini = -50;
maxi = 250;

% define plots
loadPlotParameters;

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
scatter3(xFastsum(:,1),xFastsum(:,2),pAna,5,pAna,'filled');
axis tight;
view(0, 90);
hold on;
set(gca, 'CLim', [mini maxi]);
colormap('bluewhitered');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
title('Bt (analytic)');

figure(2);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
scatter3(xFastsum(:,1),xFastsum(:,2),pCalc,5,pCalc,'filled');
axis tight;
hold on;
view(0, 90);
set(gca, 'CLim', [mini maxi]);
colormap('bluewhitered');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
title('Bt (fastsum)');

figure(3);
view(0, 90);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
axis tight;
hold on;
scatter3(xFastsum(:,1),xFastsum(:,2),pError,5,pError,'filled');
set(gca, 'CLim', [0 150]);
colormap('bluewhitered');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
title('Abs. Error (ana-fast)');


figure(4)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
scatter3(xFastsum(:,1),xFastsum(:,2),B(3,:,:),5,B(3,:,:),'filled');
axis tight;
view(0, 90);
hold on;
set(gca, 'CLim', [mini maxi]);
colormap('bluewhitered');
h = colorbar();
xlabel('x [m]');
ylabel('y [m]');
ylabel(h,'B [nT]');
title('B3 (analytic)');
hold off;