function BadPixelData = LoadBadPixelData(genum)
%	BadPixelData - Loads bad pixel data for each GE detector at 1-ID-E at
%	the APS.
%
%   INPUT:
%
%   genum
%       GE number
%
%   OUTPUT:
%
%   BadPixelData
%       list of bad pixels on the GE detector
%
%   NOTE:
%       only works for GE 1-4. Returns an empty array for other GEs.

%%% CHECK IF RUNNING AT APS SITE
if ispc
    [~, w]   = dos('HOSTNAME');
    if ~(contains(w, 'kerner') || contains(w, 'riesling') || contains(w, 'sec1parkjs'))
        disp('bad pixel data does not exist at this location.')
        disp('returning empty matrix ...')
        BadPixelData    = [];
        return
    end
elseif isunix
    [~, w]  = unix('echo $HOSTNAME');
    if ~contains(w, 'xray.aps.anl.gov')
        disp('bad pixel data does not exist at this location.')
        disp('returning empty matrix ...')
        BadPixelData    = [];
        return
    end
end

if isunix
    [~, w] = unix('whoami');
end

if contains(w, 'parkjs')
    w = 1;
elseif contains(w, 's1iduser')
    w = 2;
end

if genum == 1
    if isunix && w == 2
        fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF43522-3/Full/EF43522-3Full_BadPixel.img';
    elseif isunix && w == 1
        fname   = '/home/beams/PARKJS/ge_bad_pixel_data/EF43522-3Full_BadPixel.img';
    elseif ispc
        fname   = 'V:/misc/DetectorData/EF43522-3/Full/EF43522-3Full_BadPixel.img';
    end
    
    BadPixelData    = NreadGE(fname, 1);
    BadPixelData    = find(BadPixelData ~= 0);
elseif genum == 2
    if isunix && w == 2
        fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF44064-6/Full/EF44064-6Full_BadPixel.img';
    elseif isunix && w == 1
        fname   = '/home/beams/PARKJS/ge_bad_pixel_data/EF44064-6Full_BadPixel.img';
    elseif ispc
        fname   = 'V:/misc/DetectorData/EF44064-6/Full/EF44064-6Full_BadPixel.img';
    end
    
    BadPixelData    = NreadGE(fname, 1);
    BadPixelData    = find(BadPixelData ~= 0);
elseif genum == 3
    if isunix && w == 2
        fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF43089-5/Full/EF43089-5Full_BadPixel.img';
    elseif isunix && w == 1
        fname   = '/home/beams/PARKJS/ge_bad_pixel_data/EF43089-5Full_BadPixel.img';
    elseif ispc
        fname   = 'V:/misc/DetectorData/EF43089-5/Full/EF43089-5Full_BadPixel.img';
    end
    
    BadPixelData    = NreadGE(fname, 1);
    BadPixelData    = find(BadPixelData ~= 0);
elseif genum == 4 && w == 2
    if isunix
        fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF44066-7/Full/EF44066-7Full_BadPixel.img';
    elseif isunix && w == 1
        fname   = '/home/beams/PARKJS/ge_bad_pixel_data/EF44066-7Full_BadPixel.img';
    elseif ispc
        fname   = 'V:/misc/DetectorData/EF44066-7/Full/EF44066-7Full_BadPixel.img';
    end
    
    BadPixelData    = NreadGE(fname, 1);
    BadPixelData    = find(BadPixelData ~= 0);
elseif genum == 5 && w == 2
    if isunix
        fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/1339.6/Full/1339.6Full_BadPixel_d.txt.img';
    elseif isunix && w == 1
        fname   = '/home/beams/PARKJS/ge_bad_pixel_data/1339.6Full_BadPixel_d.txt.img';
    elseif ispc
        fname   = 'V:/misc/DetectorData/1339.6/Full/1339.6Full_BadPixel_d.txt.img';
    end
    
    BadPixelData    = NreadGE(fname, 1);
    BadPixelData    = find(BadPixelData ~= 0);
else
    disp('bad pixel data does not exist. check ge number ...')
    disp('returning empty matrix ...')
    BadPixelData    = [];
end