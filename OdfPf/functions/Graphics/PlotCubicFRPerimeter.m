function PlotCubicFRPerimeter()
% PLOTCUBICFRPERIMETER - Plot perimeter of cubic FR.
%
%   PlotCubicFRPerimeter
%
%   Plot the perimeter of the fundamental region
%   for cubic symmetries.
%
f_414  =   4.1421356200e-01; f_171  =  1.7157287500e-01;
f_414m =  -4.1421356200e-01; f_171m = -1.7157287500e-01;
%
rescale = 1.001;
linewid = 1.0;
%
cfr_top = [ f_414m f_171  f_414; ...
	    f_414m f_171m f_414; ...
	    f_171m f_414m f_414; ...
	    f_171  f_414m f_414; ...
	    f_414  f_171m f_414; ...
	    f_414  f_171  f_414; ...
	    f_171  f_414  f_414; ...
	    f_171m f_414  f_414; ...
	    f_414m f_171  f_414; ] * rescale;
%
line(cfr_top(:,1), cfr_top(:,2), cfr_top(:,3), ...
     'Color', 'k', 'LineWidth', linewid);
line(cfr_top(:,3), cfr_top(:,1), cfr_top(:,2), ...);
     'Color', 'k', 'LineWidth', linewid);
line(cfr_top(:,2), cfr_top(:,3), cfr_top(:,1), ...);
     'Color', 'k', 'LineWidth', linewid);
%
line(cfr_top(:,1), cfr_top(:,2), -1*cfr_top(:,3), ...);
     'Color', 'k', 'LineWidth', linewid);
line(-1*cfr_top(:,3), cfr_top(:,1), cfr_top(:,2), ...);
     'Color', 'k', 'LineWidth', linewid);
line(cfr_top(:,2), -1*cfr_top(:,3), cfr_top(:,1), ...);
     'Color', 'k', 'LineWidth', linewid);
%
clear f_414 f_414m f_171 f_171m cfr_top;
