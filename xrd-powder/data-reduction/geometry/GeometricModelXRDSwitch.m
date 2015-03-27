function mapped_tth = GeometricModelXRDSwitch(instr, polimg)
% GeometricModelXRDSwitch - switches between different geometric models.
%
%   USAGE:
%
%   mapped_tth = GeometricModelXRDSwitch(instr, polimg)
%
%   INPUT:
%
%   instr
%       instrument parameters in struct variable
%
%   polimg
%       radially integrated data organized in struct variable
%
%   OUTPUT:
%
%   mapped_tth
%       corrected 2theta based on the instrument
centers     = instr.centers./1000;
distance    = instr.distance;
gammaX      = instr.gammaX;
gammaY      = instr.gammaY;
detpars     = instr.detpars;
dettype     = instr.dettype;

radius      = Pixel2mm(polimg.radius', instr.pixelsize);
azimuth     = polimg.azimuth;

switch lower(dettype)
    case '0'
        disp('GeometricModelXRD0')
        mapped_tth  = GeometricModelXRD0(...
            centers, ...
            distance, ...
            gammaY, gammaX, ...
            radius, azimuth, detpars)';
    case '1'
        disp('GeometricModelXRD1')
        mapped_tth  = GeometricModelXRD1(...
            centers, ...
            distance, ...
            gammaY, gammaX, ...
            radius, azimuth, detpars)';
    case '2'
        disp('GeometricModelXRD2')
        mapped_tth  = GeometricModelXRD2(...
            centers, ...
            distance, ...
            gammaY, gammaX, ...
            radius, azimuth, detpars)';
    case '2a'
        disp('GeometricModelXRD2a')
        mapped_tth  = GeometricModelXRD2a(...
            centers, ...
            distance, ...
            gammaY, gammaX, ...
            radius, azimuth, detpars)';
    case '2b'
        disp('GeometricModelXRD2b')
        mapped_tth  = GeometricModelXRD2b(...
            centers, ...
            distance, ...
            gammaY, gammaX, ...
            radius, azimuth, detpars)';
    otherwise
        disp('Unknown method!!')
end
