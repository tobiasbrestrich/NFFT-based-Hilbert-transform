%% Load this before plotting to get optimal and uniform figures for papers

%parameters:
width = 5;     % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 12;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultAxesLineWidth',alw);   % set the default AxesLineWidth width to alw
set(0,'defaultAxesFontSize',fsz); % set the default AxesFontsize to fsz
set(0,'DefaultLegendFontSizeMode','manual');
set(0,'defaultLegendFontSize',fsz); % set the default LegendFontsize to fsz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);