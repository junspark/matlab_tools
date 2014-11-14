function Ilist=BuildMeshPolarXRD(R, V, qrule)
% numerical integration of image data

% save('BuildMeshPolarXRD_input.mat')
% return
% clear all
% close all
% clc
% 
% load('BuildMeshPolarXRD_input.mat')
% return
% qpts	= [...
%     0.44594849091597 0.44594849091597; ...
%     0.44594849091597 0.10810301816807; ...
%     0.10810301816807 0.44594849091597; ...
%     0.09157621350977 0.09157621350977; ...
%     0.09157621350977 0.81684757298046; ...
%     0.81684757298046 0.09157621350977; ...
%     ];
% qwts	= [ ...
%     0.22338158967801; ...
%     0.22338158967801; ...
%     0.22338158967801; ...
%     0.10995174365532; ...
%     0.10995174365532; ...
%     0.10995174365532; ...
%     ];

numRho  = size(V,1) - 1;
numEta  = size(V,2) - 1;

Ilist   = zeros(numRho,1);

for jj = 1:1:numRho
    vjj	= 0;
    ajj	= 0;
    for ii = 1:1:numEta
        vji     = V(jj,ii);     %%% PT1
        rji     = R(jj,ii);
        
        vji1    = V(jj,ii+1);   %%% PT2
        rji1    = R(jj,ii+1);
        
        vj1i    = V(jj+1,ii);   %%% PT3
        rj1i    = R(jj+1,ii);
        
        vj1i1   = V(jj+1,ii+1); %%% PT4
        rj1i1   = R(jj+1,ii+1);
        
        % plot(vji, rji, 'ko')
        % hold on
        % plot(vji1, rji1, 'ro')
        % plot(vj1i, rj1i, 'go')
        % plot(vj1i1, rj1i1, 'bo')
        m   = 0.5.*[ ...
            vji rji 1; ...
            vj1i rj1i 1; ...
            vji1 rji1 1; ...
            ];
        jac1    = abs(det(m));
        
        m   = 0.5.*[ ...
            vji1 rji1 1; ...
            vj1i rj1i 1; ...
            vj1i1 rj1i1 1; ...
            ];
        jac2    = abs(det(m));
        for kk = 1:1:size(qrule.pts,2)
            %%% TRIANGLE 1 (1-3-2)
            v   = vji*qrule.pts(1,kk) + vj1i*qrule.pts(2,kk) + vji1*qrule.pts(3,kk);
            r   = rji*qrule.pts(1,kk) + rj1i*qrule.pts(2,kk) + rji1*qrule.pts(3,kk);
            vjj	= vjj + v*r*qrule.wts(kk)*jac1;
            ajj	= ajj + r*qrule.wts(kk)*jac1;
            
            %%% TRIANGLE 2 (2-3-4)
            v   = vji*qrule.pts(1,kk) + vj1i1*qrule.pts(2,kk) + vj1i*qrule.pts(3,kk);
            r   = rji*qrule.pts(1,kk) + rj1i1*qrule.pts(2,kk) + rj1i*qrule.pts(3,kk);
            vjj	= vjj + v*r*qrule.wts(kk)*jac2;
            ajj = ajj + r*qrule.wts(kk)*jac2;
        end
%         return
    end
    
    Ilist(jj)   = vjj/ajj;
end