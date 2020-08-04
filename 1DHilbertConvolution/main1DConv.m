%% Computes the Gravity-Field (acceleration) of a infinite horizontal cylinder (in
%% y-direction) on a non-uniform surface Profile in x-direction.

%Add libraries and clear workspace:
clear all;
close all;
addpath('..\Util');

%%%%%%%%%%%
%% Model %%
%%%%%%%%%%%
R = 200;
N = 2001;
x = sort(randsphere(N,1,R,10));

%cylinder radius in [m]
a=5;
%depth of cylinder axix in [m] in y-direction
t=10;
%densitiy difference in [kg/m^3]
drho=500;
%center of cylinder in [m]
r0=[0,0,t];
ny=length(x);
%attraction in [mGal]
g=zeros(ny,3);

%iterate over profile
for k= 1:ny
    g(k,:)=attrHorCylinder([0,x(k),0],r0,a,drho);
end


%%%%%%%%%%%%%%%%%%%%%%
%% Frequency domain %%
%%%%%%%%%%%%%%%%%%%%%%

yyhat = linspace(-200,200,ny)';
g3hat = interp1(x,g(:,3),yyhat,'spline');
dy=uniquetol(diff(yyhat));

ky=[0:floor((ny-1)/2),-ceil((ny-1)/2):-1]'; 
ky=ky/(ny/2)*pi/dy;
hi=1i*sign(ky);
sz=fft(g3hat);
gHil1hat=real(ifft(hi.*sz));
gHilFreq = interp1(yyhat,gHil1hat,x,'spline');

%%%%%%%%%%%%%%%%%%
%% Space domain %%
%%%%%%%%%%%%%%%%%%

dx = getWeighting(x);

%Hitrafo without singularity (changed varaibles approach)
gHilCh = hiTrafo1D1NonUniform(g(:,3),x,x,dx);

%Hitrafo with near field approximation
kernel = 1./(x);
eps_I = 5;
n = round(N*eps_I/R);
gHil = zeros(ny,1);
%build approximated kernel
idx = find(abs(x-r0(2))<=eps_I);
excl = find(kernel == inf | kernel == NaN);
[A,B,YFIT] = Fseries(x,kernel,24,true,'sine');
apprKernel = kernel;
apprKernel(idx) = YFIT(idx);
gHilNF = conv(g(:,3).*dx,apprKernel,'same')./-pi;

%%%%%%%%%%%%%
%% Figures %%
%%%%%%%%%%%%
loadPlotParameters;

figure(1);clf;
hold on
plot(x, kernel./pi, 'g-','LineWidth',lw,'MarkerSize',msz);
plot(x, apprKernel./pi, 'r-','LineWidth',lw,'MarkerSize',msz);
ylabel('$\mathcal{K}(x)$', 'Interpreter','latex');
xlabel('x');
legend('$\mathcal{K}(x)$', '$\mathcal{K}_{RF}(x)$', 'Interpreter','latex', 'Orientation', 'vertical');
hold off
grid on;
set(gca,'XTick',-200:50:200); 
ylim([-0.1 0.1]);
set(gca,'YTick',-0.16:0.02:0.16); 

figure(2);clf;
hold on
plot(x, g(:,2), 'k-','LineWidth',lw,'MarkerSize',msz);
plot(x, gHilFreq, 'b-','LineWidth',lw,'MarkerSize',msz);
plot(x, gHilCh, 'g-','LineWidth',lw,'MarkerSize',msz);
plot(x,gHilNF,'m-','LineWidth',lw,'MarkerSize',msz);
ylabel('g [mGal]');
set(gca,'XTick',-200:50:200); 
ylim([-0.06 0.06]);
set(gca,'YTick',-0.06:0.01:0.06);
xlabel('x [m]');
legend('g_1', 'g_1 (freq)', 'g_1 (changed var.)','g_1 (regularized)');
hold off
grid on;