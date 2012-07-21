function RP = RPError(epsilon, refFun, fun)
% RP = RPError(epsilon, refFun, fun)

refFun = refFun(:);
fun = fun(:);

len = length(fun);

if length(refFun) ~= len
  error('funtions do not match')
end

index = find(refFun > epsilon);
theta = zeros(len, 1);
theta(index) = 1;

rij = abs(refFun(index) - fun(index))./refFun(index);

%RPij = 100*theta.*rij;
RPij = 100*rij;

RP = sum(RPij)/sum(theta);
