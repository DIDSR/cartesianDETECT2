%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE: new_ab_fit.m                         %
% AUTHORS: Diksha Sharma (US FDA, MD, USA)   %
%          Christina Sze (RMD Inc., MA, USA) %
%          Aldo Badano (US FDA, MD, USA)     %
% DATE: October 31, 2013                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% DISCLAIMER %
%%%%%%%%%%%%%%

% This software and documentation (the "Software") were developed at the
% Food and Drug Administration (FDA) by employees of the Federal Government
% in the course of their official duties. Pursuant to Title 17, Section 105
% of the United States Code, this work is not subject to copyright
% protection and is in the public domain. Permission is hereby granted,
% free of charge, to any person obtaining a copy of the Software, to deal
% in the Software without restriction, including without limitation the
% rights to use, copy, modify, merge, publish, distribute, sublicense, or
% sell copies of the Software or derivatives, and to permit persons to whom
% the Software is furnished to do so. FDA assumes no responsibility
% whatsoever for use by other parties of the Software, its source code,
% documentation or compiled executables, and makes no guarantees, expressed
% or implied, about its quality, reliability, or any other characteristic.
% Further, use of this code in no way implies endorsement by the FDA or
% confers any advantage in regulatory decisions. Although this software can
% be redistributed and/or modified freely, we ask that any derivative works
% bear some notice that they are derived from it, and any modified versions
% bear some notice that they have been modified.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain parameters a, b at various DOIs. Fit a polynomial to the curves
% and obtain 'm' and 'n' coefficients.

clear all


% use training data set to obtain a,b,c coefficients at diff DOI
ctr = 1;    
for i=0:1:9,
x=load(sprintf('radialbin%d.out',i));
maxVal = max(x(:,2)); % normalize
x(:,2) = x(:,2)/maxVal;
xdata=x(1:21,1);
ydata=x(1:21,2);
[c,g]=GaussianFit_ab(xdata,ydata);% two-gaussian distribution fit
coeff = coeffvalues(c);
cell = struct2cell(g);
A(ctr,1) = coeff(1); % height of wider gaussian - a
B(ctr,1) = coeff(2); % width of wider gaussian - b
C(ctr,1) = coeff(3); % width of narrower gaussian - c
rmse(ctr,1) = cell(5);
ctr = ctr +1;
end

% Fit a polynomial on a,b at diff depths and obtain M, N coefficients
mytemp = load('random_xyz.out');
xplot = mytemp(1:10,3)+500; % height measured from sensor plane
[linec,lineg] = FitAB_coeffM0M1M2(xplot,A); % quadratic fit for coefficient a
newcoeff = coeffvalues(linec);
cell = struct2cell(lineg);
M2 = newcoeff(1); % m2
M1 = newcoeff(2); % m1
M0 = newcoeff(3); % m0

[lineca,linega] = FitAB_coeffM0M1M2(xplot,B); % quadratic fit for coefficient b
newcoeff2 = coeffvalues(lineca);
cell2 = struct2cell(linega);
N2 = newcoeff2(1); % n2
N1 = newcoeff2(2); % n1
N0 = newcoeff2(3); % n0

hh=figure
axes1 = axes('Parent',hh,'XTickLabel',{'','','','','',''},...
    'Position',[0.366032210834553 0.583837209302326 0.321376281112738 0.341162790697675],...
    'FontSize',18);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 1000]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 1]);
box(axes1,'on');
hold(axes1,'all');

subplot(2,1,1) % plots the fit for coefficients
x1=floor(min(xplot)):01:floor(max(xplot)); %% modified
f = M2*x1.^2 + M1*x1 + M0;
plot((xplot),(A),hh,axes1,'MarkerSize',10,'Marker','*','LineStyle','none');
hold on
plot(x1,f,'LineWidth',2,'Color',[1 0 0]);
 set(gca,'XDir','reverse') % the xplot will plot the height on x-axis. reverse to obtain DOI.
axis([0 1000 0 1])
ylabel('a','FontSize',18)

subplot(2,1,2)

f = N2*x1.^2 + N1*x1 + N0;
% plot the height
plot((xplot),(B),hh,axes1,'MarkerSize',10,'Marker','*','LineStyle','none');
hold on
plot(x1,f,'LineWidth',2,'Color',[1 0 0]);
 set(gca,'XDir','reverse')
axis([0 1000 0 10])
ylabel('b','FontSize',18)
xlabel('DOI (\mum)','FontSize',18);
