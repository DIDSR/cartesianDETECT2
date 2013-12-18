%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE: new_bin_and_fit.m                    %
% AUTHOR: Diksha Sharma (US FDA, MD, USA)    % 
%         Aldo Badano (US FDA, MD, USA)      %
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

% Plots radial response and fits a Gaussian function on it

x = load('radialbin100.out'); % sample radial response array
xdata=x(1:20,1);
ydata=x(1:20,2);
ydata=ydata/max(ydata);

[c,g]=GaussianFit_ab(xdata,ydata); %Gaussian fit
figure
plot(xdata,ydata,'-b','linewidth',2); %plot the radial bins
coef = coeffvalues(c);
a=coef(1);
b=coef(2);
c=coef(3);
xi = 0:.005:20;
y= a.*exp(-(xi.^2/b.^2) ) +(1-a).*exp(-(xi.^2/c.^2));
hold on; %plot the fit
plot(xi,y,'--r','linewidth',2);