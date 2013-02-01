function mata = MatArrayOfSixVectors(vec6)
% MatArrayOfSixVectors - Array of 3x3 matrices from array of 6-vectors.\
% \'a0\'a0\
% \'a0\'a0mata = MatArrayOfSixVectors(vec6)\
%\
% \'a0\'a0vec6 is a 6 x n array of doubles,\
% \'a0\'a0\'a0\'a0\'a0\'a0\'a0representing symmetrics matrices; components\
% are 11, 21, 31, 22, 32, 33\
indices = [1 2 3 2 4 5 3 5 6];
%\
mata = reshape(vec6(indices, :), [3 3 size(vec6, 2)]);