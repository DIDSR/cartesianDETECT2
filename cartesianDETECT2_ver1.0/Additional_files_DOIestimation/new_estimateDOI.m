%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE: new_estimateDOI.m                    %
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

% Estimate DOI based on fit coefficients 'm' and 'n'
% Before executing this program, run new_ab_fit.m and do not clear the history.

myz = load('subset_random_xyz.out'); % load actual DOI data

ctr = 1;
    newctr = 0; % for correct naming the file radialbin*.out
    
for i=0:1:9,% simulations made at interval of 100 from height 30 to 97
    
    
    x=load(sprintf('radialbin%d.out',newctr)); % test simulations named in console
    newctr = newctr+1; % for naming the file radialbin*.out
        
    xdata=x(1:21,1);
    ydata=x(1:21,2);
    ydata=x(1:21,2)/max(ydata);
    
    [c,g]=GaussianFit_ab(xdata,ydata); % fits new data - obtain a1 and b1
    coeff = coeffvalues(c);
    cell = struct2cell(g);
    A1(ctr) = coeff(1); % new a
    B1(ctr) = coeff(2); % new b
    
    DOI_M_1 = (-M1 - sqrt( (M1*M1) - 4*M2*(M0-A1(ctr)) ) )/(2*M2); % using 'a1'
    DOI_M_2 = (-M1 + sqrt( (M1*M1) - 4*M2*(M0-A1(ctr)) ) )/(2*M2);
    
    DOI_N_1 = (-N1 - sqrt( (N1*N1) - 4*N2*(N0-B1(ctr)) ) )/(2*N2); % using 'b1'
    DOI_N_2 = (-N1 + sqrt( (N1*N1) - 4*N2*(N0-B1(ctr)) ) )/(2*N2);
    
    my_doiM1(ctr) = real(DOI_M_1);
    my_doiM2(ctr) = real(DOI_M_2);
    my_Z(ctr) = myz(i+1,3); % saves the actual height
    
    my_doiN1(ctr) = real(DOI_N_1);% take the real value from the height estimated from fit for b
    my_doiN2(ctr) = real(DOI_N_2);
    
    ctr = ctr +1;
end

DOI(1:ctr-1,1)=0;
DOI(1:ctr-1,2)=0;
DOI(1:ctr-1,3)=0;
% finds the height common to M and N by comparing the estimated heights from a fit and b fit. finds the smallest
% difference betweend the two and takes the average between the two.
for ctr2 = 1:1:(ctr-1),
    diff = abs(my_doiM1(ctr2) - my_doiN1(ctr2)); % starts by subtracting the first estimated height from each estimate
    if abs(my_doiM1(ctr2) - my_doiN1(ctr2)) <= diff,
        DOI(ctr2,1) = (my_doiM1(ctr2) + my_doiN1(ctr2))/2;
        diff = abs(my_doiM1(ctr2) - my_doiN1(ctr2));
    end
    if abs(my_doiM1(ctr2) - my_doiN2(ctr2)) < diff,
        DOI(ctr2,1) = (my_doiM1(ctr2) + my_doiN2(ctr2))/2;
        diff = abs(my_doiM1(ctr2) - my_doiN2(ctr2));
    end
    if abs(my_doiM2(ctr2) - my_doiN1(ctr2)) < diff,
        DOI(ctr2,1) = (my_doiM2(ctr2) + my_doiN1(ctr2))/2;
        diff = abs(my_doiM2(ctr2) - my_doiN1(ctr2));
    end
    if abs(my_doiM2(ctr2) - my_doiN2(ctr2)) < diff,
        DOI(ctr2,1) = (my_doiM2(ctr2) + my_doiN2(ctr2))/2;
        diff = abs(my_doiM2(ctr2) - my_doiN2(ctr2));
        
    end
    
    %DOI(ctr2,2) = ( DOI(ctr2,1) - (my_Z(ctr2)+500) ) / (my_Z(ctr2)+500);% takes the average of the value with the lowest difference
    DOI(ctr2,2) = ( DOI(ctr2,1) - (my_Z(ctr2)+500) );
    DOI(ctr2,3) = (my_Z(ctr2)+500); % the actual height value
    
    % ctr2 = ctr2 + 1;
    
end


hh1=figure
axes1 = axes('Parent',hh1,...
    'Position',[0.34375 0.459818880943431 0.283854166666667 0.24586306317724],...
    'FontSize',16);
%subplot(2,1,1)
plot(DOI(:,3), DOI(:,1),'MarkerSize',8,'Marker','o','LineWidth',1,'LineStyle','none');
hold on
plot(DOI(:,3),DOI(:,3),'r-','LineWidth',2);
xlabel('Actual DOI (\mum)','FontSize',18);
ylabel('Estimated DOI (\mum)','FontSize',18);

axis([0 1000 0 1000])

axes2 = axes('Parent',hh1,...
    'XTickLabel',{'1000','800','600','400','200','0'},...
    'XDir','reverse',...
    'Position',[0.34375 0.136239782016349 0.283854166666667 0.243869209809264],...
    'FontSize',16);
%subplot(2,1,2)
plot(1000-DOI(:,3),(DOI(:,2)),'MarkerSize',8,'Marker','*','LineStyle','none');
ylabel('Error (\mum)','FontSize',18)
xlabel('DOI (\mum)','FontSize',18);


%%%%%%%

