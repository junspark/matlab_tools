clear all
close all
clc

return
pname   = 'V:\stebner_dec13\hexrd_analysis\B2.LatticeConstant\corrected';
fstem   = 'NDC5';

imsum   = zeros(2048,2048);
for i = 66:1:105
    fname   = [fstem, '_', sprintf('%05d', i), '.ge2.sum'];
    pfname  = fullfile(pname, fname);
    
    fid = fopen(pfname, 'r');
    imdata  = fread(fid,[2048 2048],'*float');
    fclose(fid);
    
    imagesc(log(abs(imdata)))
    axis equal
    pause(1)
    
    imsum   = imsum + imdata;
end

pname   = 'V:\stebner_dec13\hexrd_analysis\B2.LatticeConstant\corrected';
fname   = 'NDC5_synthetic_powder.sum';
pfname  = fullfile(pname, fname);

fid = fopen(pfname, 'w');
fwrite(fid, imsum, 'float32');
fclose(fid);