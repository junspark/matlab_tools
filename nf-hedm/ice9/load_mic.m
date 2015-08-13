%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% load mic file - just like XDM Toolkit
%
%
% sidewidth is the width of one side of the triangle
%
%
% File Format:
% Col 1-3 x, y, z
% Col 4   1 = triangle pointing up, 2 = triangle pointing down
% Col 5 - Generation number
% Col 7-9 orientation
% Col 10  Confidence
%
%
function [output, sidewidth]= load_mic(filename, numCols)



snp = textread( filename );
sidewidth = snp(1, 1);
output = snp(2:end, :);

% fd = fopen(filename);
% remove first line
% sidewidth = fgetl(fd);
% sidewidth = str2double(sidewidth);
% if(numCols == 9)
%     output = fscanf(fd, '%g %g %g %d %d %d %g %g %g', [numCols, inf]);
% elseif(numCols == 10)
%     output = fscanf(fd, '%g %g %g %d %d %d %g %g %g %g', [numCols, inf]);
% elseif(numCols == 10)
%     output = fscanf(fd, '%g %g %g %d %d %d %g %g %g %g', [numCols, inf]);
% end
% output = output';

% fclose(fd);
end