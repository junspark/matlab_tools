function pfname = GenerateGEpfname(ImagePars)
% GenerateGEpfname - generate file names for GE image files
%
%   USAGE:
%
%   pfname = PlotSPF(ImagePars)
%
%   INPUT:
%
%   PFpts
%       pole figure points on the sphere (scattering vectors) (n x 3)
%
%   Data
%       scalar values corresponding to PFpts (n x 1)
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'Title'            title of the figure (default is 'SPF')
%   'ViewAngle'        angle to view the pole figure from
%                      default is [-223 26]
%   'DataRange'        range of the Data to determine the colors
%                      default is [-1 1]
%   'xaxis'            label for the x-axis
%                      default is 'x'
%   'yaxis'            label for the y-axis
%                      default is 'y'
%   'zaxis'            label for the z-axis
%                      default is 'z'
%   'ShowSurf'         on|{off}
%                      to display the PF surface
%
%   OUTPUT:  none
%
%   NOTE:
%           This function is similar to PlotPoints
%


numimages   = length(ImagePars.fnumber);
for i = 1:1:numimages
    fname   = sprintf([ImagePars.fbase, '%05d.%s', ], ImagePars.fnumber(i), ImagePars.fext);
    pfname{i,1} = fullfile(ImagePars.pname, fname);
end
