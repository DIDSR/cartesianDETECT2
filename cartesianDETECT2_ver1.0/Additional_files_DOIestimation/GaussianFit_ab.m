%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE: GaussianFit_ab.m                     %
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

function [cf_,gof] =GaussianFit_ab(xdata,ydata)
%GaussianFit_ab    Create plot of datasets and fits
%   GaussianFit_ab(XDATA,YDATA)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

 
% Data from dataset "ydata vs. xdata":
%    X = xdata:
%    Y = ydata:
%    Unweighted
%
% This function was automatically generated on 13-Jun-2011 03:50:27

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[606 244 674 479]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;

 
% --- Plot data originally in dataset "ydata vs. xdata"
xdata = xdata(:);
ydata = ydata(:);
h_ = line(xdata,ydata,'Parent',ax_,'Color',[0.333333 0 0.666667],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(xdata));
xlim_(2) = max(xlim_(2),max(xdata));
legh_(end+1) = h_;
legt_{end+1} = 'ydata vs. xdata';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-0.19000000000000000222, 19.190000000000001279]);
end


% --- Create fit "fit 1"
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[1 Inf Inf]);
ok_ = isfinite(xdata) & isfinite(ydata);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.15926779736240237462 0.76440398937158793302 0.060387460114136981559 ];
set(fo_,'Startpoint',st_);
ft_ = fittype('a*exp(-(x^2 / b^2))+(1-a)*exp(-(x^2 / c^2))',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a', 'b', 'c'});

% Fit this model using new data
cf_ = fit(xdata(ok_),ydata(ok_),ft_,fo_);

% Or use coefficients from the original fit:
if 0
   cv_ = { 0.4177916633643498634, 4.4506058366462584175, 0.060387460114136981559};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
     'LineStyle','-', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';
[cf_,gof] = fit(xdata(ok_),ydata(ok_),ft_,fo_);
% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
