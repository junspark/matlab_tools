function [mapping_table] = hedm_TrackGrains(grains0, grains1, varargin)
% hedm_TrackGrains
%   generates a mapping table between two FF-HEDM maps. 
%
%   mapping_table = hedm_TrackGrains(grains0, grains1)
%
%   INPUT:
%
%   grains0
%       grain map at the reference state. the map is read in using
%       parseGrainData or parseGrainData_OneLayer function.
%
%   grains1
%       grain map at the new state. the map is read in using
%       parseGrainData or parseGrainData_OneLayer function.
%
%   These arguments can be followed by a list of
%   parameter/value pairs. Options are:
%
%   'ang_res'               the angular resolution to be used for 
%                           weighting (deg). the default is 0.25 deg.
%   'pos_res'               the position resolution to be used for 
%                           weighting (micron). the default is 100 micron.
%   'save_mapping_table'    save mapping table to a file. the default is
%                           false.
%
%   OUTPUT:
%   mapping_table           a table containing mapping between grains0 and
%                           grains1
%
%   NOTE:
%   1. mapping is computed based on the table view of grains0 and grains1. 

% default options
optcell = {...
    'ang_res', 0.25, ...
    'pos_res', 100, ...
    'save_mapping_table', false, ...
    'mapping_table_fname', 'mapping_table.mat', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

numgrains0  = length(grains0);
numgrains1  = length(grains1);

miso_table  = zeros(numgrains0, numgrains1);
dist_table  = zeros(numgrains0, numgrains1);

grains0_COM   = [grains0(:).COM];
grains0_quat  = [grains0(:).quat];

grains1_COM   = [grains1(:).COM];   % COMii    = [grains1(:).COM];
grains1_quat  = [grains1(:).quat];  % quatii   = [grains1(:).quat];

% parpool(20)
for iii = 1:1:numgrains0
    disp(sprintf('table entries for grain number %d', iii));
    
    dist_table(iii,:)    = sqrt( ...
        (grains0_COM(1, iii) - grains1_COM(1,:)).^2 + ...
        (grains0_COM(2, iii) - grains1_COM(2,:)).^2 + ...
        (grains0_COM(3, iii) - grains1_COM(3,:)).^2 ...
        );
    miso_table(iii,:)    = rad2deg(Misorientation(grains0_quat(:,iii), grains1_quat, CubSymmetries));
end
% delete(gcp)

dist_table_master   = dist_table;
miso_table_master   = miso_table;

%%% MAPPING TABLE FORMAT
% idx_S0, idx_S1, grain_id_S0, grain_id_S1, dist, miso, reason
% mapping_table   = zeros(numgrains0, 7);
% idx_S0, idx_S1, grain_id_S0, grain_id_S1, dist, miso, reason, x0, y0,
% z0, x1, y1, z1, q01, q02, q03, q04, q11, q12, q13, q14, vol0, vol1
mapping_table   = zeros(numgrains0, 23);

tic
for i = 1:1:numgrains0
    [~, idx_dist]               = min(dist_table(:), [], 'omitnan');
    [row_S0_dist, col_S1_dist]  = ind2sub([numgrains0, numgrains1], idx_dist);
    
    [~, idx_miso]               = min(miso_table(:), [], 'omitnan');
    [row_S0_miso, col_S1_miso]  = ind2sub([numgrains0, numgrains1], idx_miso);
    
    grain_id_S0_dist    = grains0(row_S0_dist).GrainID;
    grain_id_S1_dist    = grains1(col_S1_dist).GrainID;
    grain_id_S0_miso    = grains0(row_S0_miso).GrainID;
    grain_id_S1_miso    = grains1(col_S1_miso).GrainID;
    
    if idx_dist == idx_miso
        disp(sprintf('++++++++++++++++++++++++++++++++++++'))
        disp(sprintf('found : %d in S0 == %d in S1 to match in COM and orientation', row_S0_dist, col_S1_dist))
        disp(sprintf('found : S0 grainID %d == S0 grainID %d ', grain_id_S0_dist, grain_id_S1_dist))
        disp(sprintf('min distance difference    (um) : %f', dist_table(idx_dist)))
        disp(sprintf('misorientation difference (deg) : %f', miso_table(idx_miso)))
        disp(sprintf('++++++++++++++++++++++++++++++++++++'))
        
        mapping_table(i, 1) = row_S0_dist;
        mapping_table(i, 2) = col_S1_dist;
        mapping_table(i, 3) = grain_id_S0_dist;
        mapping_table(i, 4) = grain_id_S1_dist;
        mapping_table(i, 5) = dist_table(idx_dist);
        mapping_table(i, 6) = miso_table(idx_dist);
        mapping_table(i, 7) = 0;
        mapping_table(i, 8:10)  = grains0(row_S0_dist).COM;
        mapping_table(i, 11:13) = grains1(col_S1_dist).COM;
        mapping_table(i, 14:17) = grains0(row_S0_dist).quat;
        mapping_table(i, 18:21) = grains1(col_S1_dist).quat;
        mapping_table(i, 22) = grains0(row_S0_dist).GrainRadius;
        mapping_table(i, 23) = grains1(col_S1_dist).GrainRadius;
        
        dist_table(row_S0_dist, :) = nan;
        dist_table(:, col_S1_dist) = nan;
        miso_table(row_S0_dist, :) = nan;
        miso_table(:, col_S1_dist) = nan;
    else
        disp(sprintf('++++++++++++++++++++++++++++++++++++'))
        disp(sprintf('candidate pairs'))
        disp(sprintf('++++++++++++++++++++++++++++++++++++'))
        disp(sprintf('%d in S0 and %d in S0 by distance', row_S0_dist, col_S1_dist))
        disp(sprintf('S0 grain %d and S1 grain %d by distance', grain_id_S0_dist, grain_id_S1_dist))
        disp(sprintf('min distance difference  (um) : %f', dist_table(idx_dist)))
        disp(sprintf('misorientation          (deg) : %f', miso_table(idx_dist)))
        disp(sprintf('++++++++++++++++++++++++++++++++++++'))
        disp(sprintf('%d in S0 and %d in S0 by misorientation', row_S0_miso, col_S1_miso))
        disp(sprintf('S0 grain %d and S1 grain %d by misorientation', grain_id_S0_miso, grain_id_S1_miso))
        disp(sprintf('min misorientation      (deg) : %f', miso_table(idx_miso)))
        disp(sprintf('distance difference      (um) : %f', dist_table(idx_miso)))
        disp(sprintf('++++++++++++++++++++++++++++++++++++'))
        
        %%% WEIGHTING FUNCTION
        f_dist  = dist_table(idx_dist)/opts.pos_res + miso_table(idx_dist)/opts.ang_res;
        f_miso  = dist_table(idx_miso)/opts.pos_res + miso_table(idx_miso)/opts.ang_res;
        
        if f_dist < f_miso
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            disp(sprintf('mapping based on min distance'));
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            mapping_table(i, 1) = row_S0_dist;
            mapping_table(i, 2) = col_S1_dist;
            mapping_table(i, 3) = grain_id_S0_dist;
            mapping_table(i, 4) = grain_id_S1_dist;
            mapping_table(i, 5) = dist_table(idx_dist);
            mapping_table(i, 6) = miso_table(idx_dist);
            mapping_table(i, 7) = 1;
            mapping_table(i, 8:10)  = grains0(row_S0_dist).COM;
            mapping_table(i, 11:13) = grains1(col_S1_dist).COM;
            mapping_table(i, 14:17) = grains0(row_S0_dist).quat;
            mapping_table(i, 18:21) = grains1(col_S1_dist).quat;
            mapping_table(i, 22) = grains0(row_S0_dist).GrainRadius;
            mapping_table(i, 23) = grains1(col_S1_dist).GrainRadius;
            
            dist_table(row_S0_dist, :) = nan;
            dist_table(:, col_S1_dist) = nan;
            miso_table(row_S0_dist, :) = nan;
            miso_table(:, col_S1_dist) = nan;
        else
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            disp(sprintf('mapping based on min misorientation'));
            disp(sprintf('++++++++++++++++++++++++++++++++++++'))
            mapping_table(i, 1) = row_S0_miso;
            mapping_table(i, 2) = col_S1_miso;
            mapping_table(i, 3) = grain_id_S0_miso;
            mapping_table(i, 4) = grain_id_S1_miso;
            mapping_table(i, 5) = dist_table(idx_miso);
            mapping_table(i, 6) = miso_table(idx_miso);
            mapping_table(i, 7) = 2;
            mapping_table(i, 8:10)  = grains0(row_S0_miso).COM;
            mapping_table(i, 11:13) = grains1(col_S1_miso).COM;
            mapping_table(i, 14:17) = grains0(row_S0_miso).quat;
            mapping_table(i, 18:21) = grains1(col_S1_miso).quat;
            mapping_table(i, 22)    = grains0(row_S0_miso).GrainRadius;
            mapping_table(i, 23)    = grains1(col_S1_miso).GrainRadius;
            
            dist_table(row_S0_miso, :) = nan;
            dist_table(:, col_S1_miso) = nan;
            miso_table(row_S0_miso, :) = nan;
            miso_table(:, col_S1_miso) = nan;
        end
    end
end
toc

if opts.save_mapping_table
    pfname_map  = opts.mapping_table_fname;
    pos_res = opts.pos_res;
    ang_res = opts.ang_res;
    
    disp(sprintf('saving mapping table to %s', pfname_map))
    
    dist_table  = dist_table_master;
    miso_table  = miso_table_master;
    
    save(pfname_map, 'mapping_table', 'dist_table', 'miso_table', 'ang_res', 'pos_res');
    
    % idx_S0, idx_S1, grain_id_S0, grain_id_S1, dist, miso, reason
    % mapping_table   = zeros(numgrains0, 7);
    % idx_S0, idx_S1, grain_id_S0, grain_id_S1, dist, miso, reason, x0, y0,
    % z0, x1, y1, z1, q01, q02, q03, q04, q11, q12, q13, q14, vol0, vol1
    
    %%% WRITE OUT CSV FILE AS WELL
end