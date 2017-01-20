function BadPixelData = LoadBadPixelData(genum)
%%% *.IMG FILES ARE DIFFERENT FOR EACH GE
if genum == 1
    fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF43522-3/Full/EF43522-3Full_BadPixel.img';
elseif genum == 2
    fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF44064-6/Full/EF44064-6Full_BadPixel.img';
elseif genum == 3
    fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF43089-5/Full/EF43089-5Full_BadPixel.img';
elseif genum == 4
    fname   = '/home/beams/S1IDUSER/mnt/s1a/misc/DetectorData/EF44066-7/Full/EF44066-7Full_BadPixel.img';
else
    error('bad pixel data does not exist. check ge number ...')
end
BadPixelData    = NreadGE(fname, 1);
BadPixelData    = find(BadPixelData ~= 0);
return