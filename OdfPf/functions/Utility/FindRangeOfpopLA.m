function indices = FindRangeOfpopLA(start, finish)
% FINDRANGEOFPOPLA - poop!
%   indices = FindRangeOfpopLA(start, finish)
  
num_pfs = length(start);
if size(finish, 2) ~= num_pfs
  error('Inputs do not match!')
end


indices = cell(1, num_pfs);

hasPole = start == 0;
hasAzim = finish == 90;

if ~sum(hasPole)
  for i = 1:num_pfs
    indices{i} = [(72*(start(i)/5 - 1) + 2):(72*finish(i)/5 + 1)];
  end
else
  if num_pfs > 1
    for i = find(~hasPole == 1)
      indices{i} = [(72*(start(i)/5 - 1) + 2):(72*finish(i)/5 + 1)];
    end
    for j = find(hasPole == 1)
        indices{j} = [1:(72*(finish(j)/5) + 1)];
    end
  else
    indices{1} = [1:(72*finish(1)/5 + 1)];
  end
end

