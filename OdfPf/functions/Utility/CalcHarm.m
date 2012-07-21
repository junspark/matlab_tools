function numh = CalcHarm(L, type)
if nargin == 1
    type = 'cubic';
end

if strcmp(type, 'cubic')
    repeatUnit = [1 0 0 0 1 0 1 0 1 1 1 0];
elseif strcmp(type, 'hexagonal')
    repeatUnit = [1 0 1 0 1 0];
elseif strcmp(type, 'orthorhombic')
    repeatUnit = [1 0];
end

len       = length(repeatUnit);

base      = fix((L + 1)/len);
remainder = rem(L + 1, len);

if remainder
    temp = repmat(repeatUnit, [1, base]) + reshape(repmat([0:base - 1], [len, 1]), [1, (base)*len]);
    temp = cat(2, temp, base + repeatUnit(1:remainder));
    plot([0:length(temp) - 1], temp, 'bo-')
    numh = sum(temp);
else
    temp = repmat(repeatUnit, [1, base]) + reshape(repmat([0:base - 1], [len, 1]), [1, (base)*len]);
    plot([0:length(temp) - 1], temp, 'bo-')
    numh = sum(temp);
end
