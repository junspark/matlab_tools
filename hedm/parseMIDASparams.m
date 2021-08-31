function midas_params = parseMIDASparams(pfname_midas_ps, pfimage)
% parseMIDASparams Parse the MIDAS parameter file

fid_midas_ps    = fopen(pfname_midas_ps, 'r');
while ~feof(fid_midas_ps)
    midas_ps        = fgetl(fid_midas_ps);
    if startsWith(midas_ps, 'Lsd ')
        idx     = strfind(midas_ps, ' ');
        if length(idx) == 1
            midas_params.Dsam	= str2double(midas_ps(idx+1:end));
        elseif length(idx) > 1
            midas_params.Dsam   = str2double(midas_ps(idx(1)+1:idx(2)));
        end
    elseif startsWith(midas_ps, 'OmegaFirstFile ')
        idx     = strfind(midas_ps, ' ');
        if length(idx) == 1
            midas_params.ome_ini    = str2double(midas_ps(idx+1:end));
        elseif length(idx) > 1
            midas_params.ome_ini    = str2double(midas_ps(idx(1)+1:idx(2)));
        end
    elseif startsWith(midas_ps, 'OmegaStep ')
        idx     = strfind(midas_ps, ' ');
        if length(idx) == 1
            midas_params.ome_step   = str2double(midas_ps(idx+1:end));
        elseif length(idx) > 1
            midas_params.ome_step   = str2double(midas_ps(idx(1)+1:idx(2)));
        end
    elseif startsWith(midas_ps, 'NrFilesPerSweep ')
        idx     = strfind(midas_ps, ' ');
        if length(idx) == 1
            midas_params.num_img_per_layer  = str2double(midas_ps(idx+1:end));
        elseif length(idx) > 1
            midas_params.num_img_per_layer  = str2double(midas_ps(idx(1)+1:idx(2)));
        end
    elseif startsWith(midas_ps, 'StartNr ')
        idx     = strfind(midas_ps, ' ');
        if length(idx) == 1
            midas_params.start_frame_num    = str2double(midas_ps(idx+1:end));
        elseif length(idx) > 1
            midas_params.start_frame_num    = str2double(midas_ps(idx(1)+1:idx(2)));
        end
    elseif startsWith(midas_ps, 'EndNr ')
        idx     = strfind(midas_ps, ' ');
        if length(idx) == 1
            midas_params.end_frame_num  = str2double(midas_ps(idx+1:end));
        elseif length(idx) > 1
            midas_params.end_frame_num  = str2double(midas_ps(idx(1)+1:idx(2)));
        end
    elseif startsWith(midas_ps, 'BC ')
        idx     = strfind(midas_ps, ' ');
        midas_params.CENX  = str2double(midas_ps(idx(1)+1:idx(2)));
        midas_params.CENY  = str2double(midas_ps(idx(2)+1:idx(3)));
    end
end
fclose(fid_midas_ps);

% if midas_params.num_img_per_layer ~= length(pfimage)
%     error('number of files per layer and files designated inconsistent')
% else
%     disp('number of files per layer and files designated consistent')
% end

midas_params.total_numframes    = CalcNumFramesGE(pfimage);
% for iii = 1:1:midas_params.num_img_per_layer
%     midas_params.total_numframes    = CalcNumFramesGE(pfimage{iii});
% end

midas_params.num_frames     = midas_params.end_frame_num - midas_params.start_frame_num + 1;
midas_params.junk_frames    = midas_params.total_numframes - midas_params.num_frames;
midas_params.ome_fin        = midas_params.ome_ini + midas_params.num_frames*midas_params.ome_step;
midas_params.ome_grid       = linspace(midas_params.ome_ini, midas_params.ome_fin, midas_params.num_frames);
midas_params.pfname         = pfname_midas_ps;