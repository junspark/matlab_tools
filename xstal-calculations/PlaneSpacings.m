function [d, th] = PlaneSpacings(latticeParms, latticeType, hkls, varargin)
%   USAGE:
%   [d] = PlaneSpacings(latticeParms, latticeType, hkls)
%   [d] = PlaneSpacings(latticeParms, latticeType, hkls, lambda)
%   [d, th] = PlaneSpacings(latticeParms, latticeType, hkls, lambda)
%
%   INPUT:
%
%   latticeParms
%       lattice parameter of a unit cell in angstroms
%
%   lattType
%       lattice type
%
%   hkls
%       list of hkls (3 x n array)
%
%   lambda
%       x-ray wavelength in angstrom
%
%   OUTPUT:
%
%   d
%       plane spacing for hkls (1 x n arraym angstrom)
%
%   th
%       bragg angles for hkls (1 x n array, deg)
%       if lambda is not provided, lambda is set to 1.5418 angstroms.
%
%   forbidden: BCC h+k+l odd
%              FCC h, k, l mix of even and odd
%              DC  h, k, l even with h+k+l ~= 4n
%              HCP h+2k = 3n with l odd
if nargin == 4
    lambda = varargin{1};
end

if strcmp(latticeType, 'cubic')
    a = latticeParms;
    for i = 1:size(hkls, 2)
        d(i)  = a/(sqrt(hkls(1, i)^2 + hkls(2, i)^2 + hkls(3, i)^2));
    end
elseif strcmp(latticeType, 'hexagonal') | strcmp(latticeType, 'trigonal')
    a = latticeParms(1);
    c = latticeParms(2);
    if size(hkls, 1) == 4
        hkls = ConvertMillerBravais(hkls, 'contract');
    end
    for i = 1:size(hkls, 2)
        oneBydsq = ...
            4/3*(hkls(1, i)^2 + hkls(1, i)*hkls(2, i) + hkls(2, i)^2)/a^2 + ...
            hkls(3, i)^2/c^2;
        d(i)  = sqrt(1/oneBydsq);
    end
elseif strcmp(latticeType, 'tetragonal')
    a = latticeParms(1);
    c = latticeParms(2);
    for i = 1:size(hkls, 2)
        oneBydsq = ...
            (hkls(1, i)^2 + hkls(2, i)^2)/a^2 + ...
            hkls(3, i)^2/c^2;
        d(i)  = sqrt(1/oneBydsq);
    end
elseif strcmp(latticeType, 'orthorhombic')
    a = latticeParms(1);
    b = latticeParms(2);
    c = latticeParms(3);
    for i = 1:size(hkls, 2)
        oneBydsq = ...
            hkls(1, i)^2/a^2 + ...
            hkls(2, i)^2/b^2 + ...
            hkls(3, i)^2/c^2;
        d(i)  = sqrt(1/oneBydsq);
    end
else
    error('Unrecognized lattice type.')
end

if nargout > 1
    if nargin == 3
        lambda = 1.5418;
    end
    th = rad2deg(asin(lambda./(2*d)));
end