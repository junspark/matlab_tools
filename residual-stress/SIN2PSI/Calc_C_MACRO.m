clear all
close all
clc

angle_pos   = [0 90]';
radial_pos  = 1:1:21;

na  = length(angle_pos);
nr  = length(radial_pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE DIRECTIONAL MODULUS
%%% e11, e22, e33, e12, e13, e23
%%% SHEARS ARE 2*HOSFORD
%%% Ti811
C   = [ ...
    141  76.9 57.9 0      0      0; ...
    76.9 141  57.9 0      0      0; ...
    59.7 57.9 163  0      0      0; ...
    0    0    0    2*32.1 0      0; ...
    0    0    0    0      2*48.7 0; ...
    0    0    0    0      0      2*48.7; ...
    ];
load wshex2x

% %%% TEST CASE - ISOTROPIC
% C   = [ ...
%     248 155 155   0   0   0; ...
%     155 248 155   0   0   0; ...
%     155 155 248   0   0   0; ...
%     0    0    0   2*124  0   0; ...
%     0    0    0   0   2*124  0; ...
%     0    0    0   0   0   2*124; ...
%     ];

% %%% TEST CASE - CUBIC
% C   = [ ...
%     522 204 204   0   0   0; ...
%     204 522 204   0   0   0; ...
%     204 204 522   0   0   0; ...
%     0    0    0   2*161  0   0; ...
%     0    0    0   0   2*161  0; ...
%     0    0    0   0   0   2*161; ...
%     ];
% load wscub
% wshex2x = wscub;

RMat    = RMatOfQuat(QuatOfRod(wshex2x.frmesh.crd));
for i = 1:1:wshex2x.frmesh.numind
    % [RRt6(:,:,i), RtR6(:,:,i)]  = VoigtCOB(RMat(:,:,i));
    T(:,:,i)    = VectorizedCOBMatrix(RMat(:,:,i));
end
for j = 1:1:nr
    for i = 1:1:na
        pname_odf   = '/home/jspark/Documents/work/prjResidualStress/Ti-RS/Atish/Stress_code-2D/Final_ODF';
        fname_odf   = ['odf_Alp', num2str(angle_pos(i)), ...
            '_pos', num2str(radial_pos(j)), '.mat'];
        pfname_odf  = fullfile(pname_odf, fname_odf);
        odf = load(pfname_odf);
        odf = odf.odf;
        
        odf = ones(wshex2x.frmesh.numind,1);    % ISOTROPIC TEXTURE (TYPICAL SIN2PSI APPROACH)
        
        for k = 1:1:wshex2x.frmesh.numind
            % Ck(:,:,k)   = RRt6(:,:,k)*C*RtR6(:,:,k);
            Ck(:,:,k)   = T(:,:,k)'*C*T(:,:,k);
        end
        
        for mm = 1:1:6
            for nn = 1:1:6
                C_MACRO{i,j}(mm,nn) = odf'*wshex2x.frmesh.l2ip*squeeze(Ck(mm,nn,:));
                C_MACRO{i,j}(mm,nn) = C_MACRO{i,j}(mm,nn)/sum(odf'*wshex2x.frmesh.l2ip);
            end
        end
        
        %%% SYMMETRIZE (THERE ARE SOME NUMERICAL ROUNDING ERROR 
        %%% C_MACRO{i,j} AS IS COMES OUT SLIGHTLY ASYMMETRIC
        C_MACRO{i,j}    = (C_MACRO{i,j} + C_MACRO{i,j}')/2;
    end
end