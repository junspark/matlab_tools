function [f] = PFI_nrmFac_func(xi, odfpf, pfData, l2ip, h1sip, beta, numh, V)

global ODF

numPfs = length(pfData.data);

sclPfData = pfData;
for i = 1:numPfs
    sclPfData.data{i} = xi(i)*pfData.data{i};
end

if nargin == 6
    [tempODF, repfs] = InvertPoleFigComp(...
        odfpf, sclPfData, ...
        [], l2ip, h1sip, ...%         'InitODF', ODF, ...
        'H1Weight', beta, ...
        'Display', 'off');
elseif nargin == 7
    [tempODF, repfs] = InvertPoleFigComp(...
        odfpf, sclPfData, ...
        [], l2ip, h1sip, ...%         'InitODF', ODF, ...
        'H1Weight', beta, ...
        'Display', 'off', ...
        'DiscHarm', numh);
elseif nargin == 8
    [tempODF, repfs] = InvertPoleFigComp(...
        odfpf, sclPfData, ...
        [], l2ip, h1sip, ...%         'InitODF', ODF, ...
        'H1Weight', beta, ...
        'Display', 'off', ...
        'DiscHarm', numh, ...
        'EigMatrix', V);
end

% Assign to global explicitly
if min(tempODF) < 0 | MeanValue(tempODF, l2ip) < (1 - 1e-4) | MeanValue(tempODF, l2ip) > (1 + 1e-4)
    error('What the hell is going on!!!')
end
ODF = tempODF;

f = 0;
for i = 1:numPfs
    %% f = f + (repfs.data{i} - sclPfData.data{i})'*(repfs.data{i} - sclPfData.data{i});
    f = f + ...
        repfs.data{i}'*repfs.data{i} ...
        - 2*sclPfData.data{i}'*repfs.data{i} ...
        + sclPfData.data{i}'*sclPfData.data{i};
end
% disp(xi')
% disp(f)