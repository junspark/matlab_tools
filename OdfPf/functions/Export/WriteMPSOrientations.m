function header = WriteMPSOrientations(filename,orient, hinfo)
% WriteMPSOrientations - Write orientations in MPS format.
%   
%   USAGE:
%
%   header = WriteMPSOrientations(filename,orient, hinfo)
%
%   INPUT:
%
%   filename is a string, 
%            the name of the output file
%   orient   is m x n, 
%            the array of orientations to write
%   hinfo    is a cell array,
%            it contains header information; options are:
% 
%   'Weights',          array
%
%   'Hardnesses',       array
%
%   'Parameterization', string  [required]
%
%   'Angles',           string  ('degrees'|'radians')
%
%   OUTPUT:
%
%   header is a cell array of strings,
%          it contains the MPS (material point simulator) 
%          header information; the header and data are
%          written to the requested file
%
[m, n] = size(orient);
%
strnum = '%number-of-orientations    ';
strwts = '%read-weights              ';
strhrd = '%read-hardnesses           ';
strang = '%angle-units               ';
strprm = '%parameterization          ';
strend = '%end-header';
%
valnum = int2str(n);
valwts = 'false';
valhrd = 'false';
valang = 'radians';
valprm = 'NONE SPECIFIED';
%
nopt = length(hinfo)/2;
hopt = reshape(hinfo, [2 nopt]);
%
for opt=1:nopt
  %
  str = hopt{1, opt};
  val = hopt{2, opt};
  %
  if (strcmp(str, 'Weights'))
    valwts = 'true';
    wts    = val;
    continue
  end
  %
  if (strcmp(str, 'Hardnesses'))
    valhrd = 'true';
    hrd = val;
    continue
  end
  %
  if (strcmp(str, 'Parameterization'))
    valprm = val;
    continue
  end
  %
  if (strcmp(str, 'Angles'))
    valang = val;
    continue
  end
  %
  %  Pass unrecognized options to WriteDataFile. (to be done)
  %
end
%
if (strcmp(valwts, 'true'))
  orient = [orient; wts];
end
%
if (strcmp(valhrd, 'true'))
  orient = [orient; hrd];
end
%
header = {...
    [strnum valnum], [strwts valwts], [strhrd valhrd], ...
    [strang valang], [strprm valprm], strend           ...
         };
%
WriteDataFile(filename, header, orient');
