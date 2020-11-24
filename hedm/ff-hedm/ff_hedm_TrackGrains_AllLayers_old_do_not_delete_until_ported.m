clear all
close all
clc

%%% CREATE THE MASTER FILE
% LoadNumber  = [0, 1, 2, 4, 5, 7, 8];
% for i = 1:1:length(LoadNumber)
%     fname   = sprintf('HeatHTNS9_S%d.mat', LoadNumber(i))
%     data(i)    = load(fname);
% end
% save('HeatHTNS9.mat', 'data', 'LoadNumber');
% return

ang_res = 0.25; % deg
pos_res = 100;  % um

wsname              = 'wscub3x';      % workspace name

samplename          = 'hr1_s1';
phasename           = 'fcc';
pname_midas_root    = '/net/wolf/data/tomo1/mpe1_oct20_midas/ff/mishra_CuHEA_fcc';
pfname_xlsx         = '/net/wolf/data/tomo1/mpe1_oct20_midas/ff/mishra_CuHEA_fcc/mishra_CuHEA_fcc.xlsx';
ref_load_line_num   = 1:5;

y_offset_between_layer  = 200;

%%%%% STARTS HERE
load(wsname);
eval(['ws   = ', wsname, ';']);
clear(wsname)

num_layers_per_hedm_data_pt = length(ref_load_line_num);

table_xlsx          = readtable(pfname_xlsx, 'FileType', 'spreadsheet', 'sheet', samplename, 'ReadVariablenames', 0);
num_hedm_data_pts   = size(table_xlsx);
if rem(num_hedm_data_pts(1,1), num_layers_per_hedm_data_pt) == 0
    disp('all layers seem to be accounted for. continuing with tracking table generation.')
else
    disp('check midas outputs')
    return
end

num_hedm_data_pts   = num_hedm_data_pts(1,1)/num_layers_per_hedm_data_pt;

for iii = 1:1:num_layers_per_hedm_data_pt
    piii        = table2cell(table_xlsx(iii, 2));
    piii        = piii{1};
    
    pfname_Grains_dot_csv_at_S0{iii,1}  = fullfile(pname_midas_root, piii);
end

%%% READ GRAINS.CSV FROM S0
Grains_dot_csv_data_at_S0   = parseGrainData(pfname_Grains_dot_csv_at_S0, ws.frmesh.symmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'OffsetDirection', 'y', ...
    'OffsetValue', y_offset_between_layer);
numgrains_at_S0 = length(Grains_dot_csv_data_at_S0);
COM0            = [Grains_dot_csv_data_at_S0(:).COM];
quat0           = [Grains_dot_csv_data_at_S0(:).quat];

for iii = 1:1:num_hedm_data_pts
    for jjj = 1:1:num_layers_per_hedm_data_pt
        idx     = (iii-1)*num_layers_per_hedm_data_pt + jjj;
        pjjj    = table2cell(table_xlsx(idx, 2));
        pjjj    = pjjj{1};
        
        pfname_Grains_dot_csv_at_S1{jjj,1}  = fullfile(pname_midas_root, pjjj);
    end
    fname_map   = sprintf('mapping_%s_%s_S0_to_S%d.mat', samplename, phasename, iii);
    
    %%% READ GRAINS.CSV FROM SN
    Grains_dot_csv_data_at_S1   = parseGrainData(pfname_Grains_dot_csv_at_S1, ws.frmesh.symmetries, ...
        'CrdSystem', 'APS', ...
        'LabToSample', 0, ...
        'OffsetDirection', 'y', ...
        'OffsetValue', y_offset_between_layer);
    numgrains_at_S1 = length(Grains_dot_csv_data_at_S1);
    
    miso_table  = zeros(numgrains_at_S0, numgrains_at_S1);
    dist_table  = zeros(numgrains_at_S0, numgrains_at_S1);
    
    COM1            = [Grains_dot_csv_data_at_S1(:).COM];
    quat1           = [Grains_dot_csv_data_at_S1(:).quat];
    
    %%% MAKE MISO AND DIST TABLES
%     parpool(20)
    for jjj = 1:numgrains_at_S0
        disp(sprintf('table entries for grain number %d', jjj));
        COMjjj  = Grains_dot_csv_data_at_S0(jjj).COM;
        quatjjj = Grains_dot_csv_data_at_S0(jjj).quat;
        
        dist_table(jjj,:)    = sqrt((COMjjj(1) - COM1(1,:)).^2 + (COMjjj(2) - COM1(2,:)).^2 + (COMjjj(3) - COM1(3,:)).^2);
        miso_table(jjj,:)    = rad2deg(Misorientation(quatjjj, quat1, HexSymmetries));
    end
%     delete(gcp)
    
    dist_table_master   = dist_table;
    miso_table_master   = miso_table;
    
    dist_table_vec  = dist_table(:);
    miso_table_vec  = miso_table(:);
    [dist_table_vec_sorted, dist_table_vec_idx] = sort(dist_table_vec);
    [miso_table_vec_sorted, miso_table_vec_idx] = sort(miso_table_vec);
    
    %%% FIND THE PAIRS BETWEEN S0 AND S1
    mapping_table   = zeros(numgrains_at_S0, 7);   % idx_S0, idx_S1, grain_id_S0, grain_id_S1, dist, miso, reason
    
    tic
    for jjj = 1:1:numgrains_at_S0
        [~, idx_dist]               = min(dist_table(:), [], 'omitnan');
        [idx_S0_dist, idx_S1_dist]  = ind2sub([numgrains_at_S0, numgrains_at_S1], idx_dist);
        
        [~, idx_miso]               = min(miso_table(:), [], 'omitnan');
        [idx_S0_miso, idx_S1_miso]  = ind2sub([numgrains_at_S0, numgrains_at_S1], idx_miso);
        
        grain_id_S0_dist    = Grains_dot_csv_data_at_S0(idx_S0_dist).GrainID;
        grain_id_S1_dist    = Grains_dot_csv_data_at_S1(idx_S1_dist).GrainID;
        grain_id_S0_miso    = Grains_dot_csv_data_at_S0(idx_S0_miso).GrainID;
        grain_id_S1_miso    = Grains_dot_csv_data_at_S1(idx_S1_miso).GrainID;
        
        if idx_dist == idx_miso
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            disp(sprintf('found : %d in S0 == %d in S1 to match in COM and orientation', idx_S0_dist, idx_S1_dist))
            disp(sprintf('found : S0 grainID %d == S0 grainID %d ', grain_id_S0_dist, grain_id_S1_dist))
            disp(sprintf('min distance difference    (um) : %f', dist_table(idx_dist)))
            disp(sprintf('misorientation difference (deg) : %f', miso_table(idx_miso)))
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            
            mapping_table(jjj, 1) = idx_S0_dist;
            mapping_table(jjj, 2) = idx_S1_dist;
            mapping_table(jjj, 3) = grain_id_S0_dist;
            mapping_table(jjj, 4) = grain_id_S1_dist;
            mapping_table(jjj, 5) = dist_table(idx_dist);
            mapping_table(jjj, 6) = miso_table(idx_dist);
            mapping_table(jjj, 7) = 0;
            
            dist_table(idx_S0_dist, :) = nan;
            dist_table(:, idx_S1_dist) = nan;
            miso_table(idx_S0_dist, :) = nan;
            miso_table(:, idx_S1_dist) = nan;
        else
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            disp(sprintf('candidate pairs'))
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            disp(sprintf('%d in S0 and %d in S0 by distance', idx_S0_dist, idx_S1_dist))
            disp(sprintf('S0 grain %d and S1 grain %d by distance', grain_id_S0_dist, grain_id_S1_dist))
            disp(sprintf('min distance difference  (um) : %f', dist_table(idx_dist)))
            disp(sprintf('misorientation          (deg) : %f', miso_table(idx_dist)))
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            disp(sprintf('%d in S0 and %d in S0 by misorientation', idx_S0_miso, idx_S1_miso))
            disp(sprintf('S0 grain %d and S1 grain %d by misorientation', grain_id_S0_miso, grain_id_S1_miso))
            disp(sprintf('min misorientation      (deg) : %f', miso_table(idx_miso)))
            disp(sprintf('distance difference      (um) : %f', dist_table(idx_miso)))
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            
            f_dist  = dist_table(idx_dist)/pos_res + miso_table(idx_dist)/ang_res;
            f_miso  = dist_table(idx_miso)/pos_res + miso_table(idx_miso)/ang_res;
            
            if f_dist < f_miso
                disp(sprintf('++++++++++++++++++++++++++++++++++++'))
                disp(sprintf('mapping based on min distance'));
                disp(sprintf('++++++++++++++++++++++++++++++++++++'))
                mapping_table(jjj, 1) = idx_S0_dist;
                mapping_table(jjj, 2) = idx_S1_dist;
                mapping_table(jjj, 3) = grain_id_S0_dist;
                mapping_table(jjj, 4) = grain_id_S1_dist;
                mapping_table(jjj, 5) = dist_table(idx_dist);
                mapping_table(jjj, 6) = miso_table(idx_dist);
                mapping_table(jjj, 7) = 1;
                
                dist_table(idx_S0_dist, :) = nan;
                dist_table(:, idx_S1_dist) = nan;
                miso_table(idx_S0_dist, :) = nan;
                miso_table(:, idx_S1_dist) = nan;
            else
                disp(sprintf('++++++++++++++++++++++++++++++++++++'))
                disp(sprintf('mapping based on min misorientation'));
                disp(sprintf('++++++++++++++++++++++++++++++++++++'))
                mapping_table(jjj, 1) = idx_S0_miso;
                mapping_table(jjj, 2) = idx_S1_miso;
                mapping_table(jjj, 3) = grain_id_S0_miso;
                mapping_table(jjj, 4) = grain_id_S1_miso;
                mapping_table(jjj, 5) = dist_table(idx_miso);
                mapping_table(jjj, 6) = miso_table(idx_miso);
                mapping_table(jjj, 7) = 2;
                
                dist_table(idx_S0_miso, :) = nan;
                dist_table(:, idx_S1_miso) = nan;
                miso_table(idx_S0_miso, :) = nan;
                miso_table(:, idx_S1_miso) = nan;
            end
        end
    end
    toc
    save(fname_map, 'mapping_table', 'dist_table', 'miso_table', 'ang_res', 'pos_res', ...
        'Grains_dot_csv_data_at_S1', 'Grains_dot_csv_data_at_S0', 'pfname_Grains_dot_csv_at_S1', 'pfname_Grains_dot_csv_at_S0');
end