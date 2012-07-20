function pfname = GenerateGEpfname(img,i,j)
nn  = '00000';

pname   = img.pname;
fbase   = img.fbase;

fnumber	= img.fnumber(i,j);
% fnumber = img.fnumber(i);
fnumber = num2str(fnumber);

nn(6-length(fnumber):end) = fnumber;

fname   = [fbase nn];

pfname  = fullfile(pname, fname);