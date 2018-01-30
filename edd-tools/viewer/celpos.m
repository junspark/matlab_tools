% calculate box position in a frame
%

function [pos] = celpos(col,row,box,edge)

if ~exist('box','var'); box = [120 30];end
if ~exist('edge','var'); edge = [10 10];end

pos = [(col-1)*(box(1)+edge(1)) (row-1)*(box(2)+edge(1))];


end
