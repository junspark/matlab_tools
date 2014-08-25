function status = GenerateSyntheticPowder(pname, fstem, imnum, fext, fout)

% pname   = 'V:\stebner_dec13\hexrd_analysis\B2.LatticeConstant\corrected';
% fstem   = 'NDC5';
% fext    = '.ge2.sum';
% fout    = 'NDC5_synthetic_powder.sum';
% imnum   = 66:1:105;

imsum   = zeros(2048,2048);
for i = 1:1:length(imnum)
    fname   = [fstem, '_', sprintf('%05d', imnum(i)), fext];
    pfname  = fullfile(pname, fname);
    
    imi = ReadSUM(pfname);
    
    imagesc(log(abs(imi)))
    axis equal
    pause(1)
    
    imsum   = imsum + imi;
end

pfname  = fullfile(pname, fout);
status  = WriteSUM(pfname, imsum);