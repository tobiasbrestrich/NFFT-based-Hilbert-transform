%% Computes the Gravity-Field (acceleration) of a infinite horizontal cylinder (in
%% y-direction) on a non-uniform surface Profile in x-direction.
%% uses the Fastsum-Method and Compares its error and time with its direct
%% copmutation (Convolution) and the same calculation in frequency domain
%% (calculated on a uniform profile, derived by interpolation)

%Add libraries and clear workspace:
addpath('..\nfft-3.5.1-mexw64-openmp\fastsum');
addpath('..\Util');

clear all;
close all;

%%%%%%%%%%%
%% Model %%
%%%%%%%%%%%
si2mgal = 1e5;
eps_B = 1/160;   % outer boundary
scaling = 100;   % profile must be on -0.25 < y < 0.25; scaling can be used to get realistical dimensions

dy = 0.0001;
x = ((-0.25+eps_B/2):dy:(0.25-eps_B/2))*scaling;

%cylinder radius in [m]
a=100/scaling;
%depth of cylinder axix in [m] in y-direction
t=100/scaling;
%densitiy difference in [kg/m^3]
drho=500;
%center of cylinder in [m]
r0=[0,0.001,t];
ny=length(x);
%attraction in [mGal]
g=zeros(ny,3);

%iterate over profile
for k= 1:ny
        g(k,:)=attrHorCylinder([0,x(k),0],r0,a,drho);
end

%%%%%%%%%%%%%
%% FASTSUM %%
%%%%%%%%%%%%%

%% Initialize parameters
d = 1;          % number of dimensions
N = length(x)-1;   % number of source knots
M = N;              % number of target knots
kernel = 'one_over_x';
c = 1;          % kernel parameter
p = 8;          % degree of smoothness of regularization
flags = 1;      % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
n = 600;        % expansion degree
eps_I = p/n*1e-10;    % inner boundary
eps_B = eps_B;   % outer boundary
m = p;
nn_oversampled=2*n; % oversampling factor for NFFT



%% random source nodes
ri1 = unique(sort(ceil(randsphere(N,1,1,1)*length(x)/2+length(x)/2)));
xFastsum = x(ri1)'/scaling;
weightedAreas = getWeighting(xFastsum);
N = length(ri1);
%coefficient
alpha = g(ri1,3).*weightedAreas;
%% random target nodes
ri2 = ri1; %target nodes are the same as the source nodes
M = length(ri2);
yFastsum = x(ri2)'/scaling;

%% Perform fastsum
plan=fastsum_init(d,kernel,c,flags,n,p,eps_I,eps_B);
fastsum_set_x(plan,xFastsum,nn_oversampled,m)
fastsum_set_alpha(plan,alpha)
fastsum_set_y(plan,yFastsum,nn_oversampled,m)

tic;
fastsum_trafo_direct(plan)   % direct computation
g1Direct = fastsum_get_f(plan);
tDirect = toc;

tic;
fastsum_trafo(plan)         % fast computation
g1Fastsum = fastsum_get_f(plan);
tFastsum = toc;
fastsum_finalize(plan)


%% post calculations
g1Fastsum = g1Fastsum/-pi;
g1Direct = g1Direct/-pi;


%%%%%%%%%%%%%%%
%% Frequency %%
%%%%%%%%%%%%%%%
tic;
NInterpol = N;
xInterpol = linspace(x(1)/scaling,x(end)/scaling,NInterpol)';
dyInterpol = min(unique(diff(xInterpol)));
g3Interpol = interp1(xFastsum,g(ri1,3),xInterpol,'linear','extrap');    

%% Field transform (gx from gz)  
%Frequency domain:
ky = [0:floor((NInterpol-1)/2),-ceil((NInterpol-1)/2):-1]'; 
ky = ky/(NInterpol/2)*pi/dyInterpol;
hi = 1i*sign(ky);
sz = fft(g3Interpol);
g1Hat = real(ifft(hi.*sz));
g1Freq = interp1(xInterpol,g1Hat,yFastsum);
tFreq = toc;

%%%%%%%%%%%%%%%%
%% Statistics %%
%%%%%%%%%%%%%%%%
%Fastsum
relErrFast = abs((g(ri2,2)-g1Fastsum)./g(ri2,2));
maxErrFast = max(abs((g(ri2,2)-g1Fastsum))./abs(g(ri2,2)));
%Direct
relErrDir = abs((g(ri2,2)-g1Direct)./g(ri2,2));
maxErrDir = max(abs((g(ri2,2)-g1Direct))./abs(g(ri2,2)));
%Interpolation
relErrInterp = abs((g(ri2,2)-g1Freq)./g(ri2,2));
maxErrInterp = max(abs((g(ri2,2)-g1Freq))./abs(g(ri2,2)));


%% plot source and target evaluations
close all;

figure(1)
hold on
p1=plot(x(ri2), g(ri2,2), 'k');
p2=plot(x(ri2), real(g1Fastsum), 'r--');
p3=plot(x(ri2), g1Freq, 'b--');
ylabel('g in [mGal]');
% yyaxis right
% p3=plot(x(ri2), rel_err_fastsum, 'r.-');
% ylabel('rel. Error');
xlabel('Profil in [m]');
legend('g1', 'g1 (fastsum)', 'g1 (frequency)');
hold off