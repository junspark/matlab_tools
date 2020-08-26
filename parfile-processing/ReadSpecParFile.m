function pardata = ReadSpecParFile(fname, varargin)
% ReadSpecParFile - read par file generated from python
%
%   INPUT:
%
%   fname
%       name of the spec generated par file name
%
%   OUTPUT:
%
%   pardata
%       spec par file data in struct arrary.
%
%   The input arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'numchs'
%       number of channels recorded in the spec generated par file (excluding the
%       date / time channel) for arbitrary format
%
%   'Version'
%       version of spec metadata file (par file).
%       example: (coratella_feb15). default value is
%           'none'. 'mpe_standard' is the new MPE metadata format for 
%           1-BM-B EDD / 6-BM-A EDD / 1-ID data.
%       this option overrides 'numchs' input.

% default options
optcell = {...
    'Version', 'none', ...
    'NumChannels', 58, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

%%% MAKE FORMAT STRING
switch lower(opts.Version)
    case 'none'
        
        disp('user specified par file format')
        numchs      = opts.NumChannels;
        
        fmtstring   = '%s %s %s %s %d %s';
        for i = 1:1:numchs
            fmtstring   = [fmtstring, ' %f'];
        end
        
        %%% READ IN DATA USING FORMAT STRING
        fid         = fopen(fname, 'r','n');
        textdata    = textscan(fid, fmtstring);
        fclose(fid);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        pardata.froot   = textdata{6};
        pardata.chs     = zeros(length(pardata.date), numchs);
        for i = 1:1:numchs
            pardata.chs(:,i)  = textdata{i+6};
        end
    case 'saxs_waxs_fmt_fastpar_v1'
        %%% STARTING internal_jun19
        fmtstring   = ['%s %s %d %s %d ' ...
            '%u %f %f %f %f %f %d %d ' ....
            '%s %s %f %f %d %s %d %d ' ...
            '%s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f' ...
            '%f %f %f %f'];
        
        %%% READ IN DATA USING FORMAT STRING
        fid     = fopen(fname, 'r','n');
        textdata  = textscan(fid, fmtstring);
        fclose(fid);
        
        whos textdata
        whos fmtstring
        %%% PARSE DATA
        %%% LINE 1
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        %%% LINE 2
        pardata.epoch_time  = textdata{6};
        pardata.integ_time  = textdata{7};
        pardata.Iring       = textdata{8};
        pardata.und_gap     = textdata{9};
        pardata.energy      = textdata{10};
        pardata.energy_cal  = textdata{11};
        pardata.foil_pos    = textdata{12};
        pardata.atten_pos   = textdata{13};
        
        %%% LINE 3
        pardata.det_type        = textdata{14};
        pardata.scan_mtr        = textdata{15};
        pardata.scan_ini        = textdata{16};
        pardata.scan_fin        = textdata{17};
        pardata.scan_nframes    = textdata{18};
        pardata.imgprefix       = textdata{19};
        pardata.imnum_ini       = textdata{20};
        pardata.imnum_fin       = textdata{21};
        
        %%% LINE 4
        pardata.det1_fname              = textdata{22};
        pardata.det1_fnum               = textdata{23};
        pardata.det1_frames_per_file    = textdata{24};
        pardata.det1_time_per_frame     = textdata{25};
        
        pardata.det2_fname              = textdata{26};
        pardata.det2_fnum               = textdata{27};
        pardata.det2_frames_per_file    = textdata{28};
        pardata.det2_time_per_frame     = textdata{29};
        
        pardata.det3_fname              = textdata{30};
        pardata.det3_fnum               = textdata{31};
        pardata.det3_frames_per_file    = textdata{32};
        pardata.det3_time_per_frame     = textdata{33};
        
        pardata.det4_fname              = textdata{34};
        pardata.det4_fnum               = textdata{35};
        pardata.det4_frames_per_file    = textdata{36};
        pardata.det4_time_per_frame     = textdata{37};
        
        pardata.det5_fname              = textdata{38};
        pardata.det5_fnum               = textdata{39};
        pardata.det5_frames_per_file    = textdata{40};
        pardata.det5_time_per_frame     = textdata{41};
        
        pardata.det6_fname              = textdata{42};
        pardata.det6_fnum               = textdata{43};
        pardata.det6_frames_per_file    = textdata{44};
        pardata.det6_time_per_frame     = textdata{45};
        
        pardata.det7_fname              = textdata{46};
        pardata.det7_fnum               = textdata{47};
        pardata.det7_frames_per_file    = textdata{48};
        pardata.det7_time_per_frame     = textdata{49};
        
        pardata.det8_fname              = textdata{50};
        pardata.det8_fnum               = textdata{51};
        pardata.det8_frames_per_file    = textdata{52};
        pardata.det8_time_per_frame     = textdata{53};
        
        pardata.det9_fname              = textdata{54};
        pardata.det9_fnum               = textdata{55};
        pardata.det9_frames_per_file    = textdata{56};
        pardata.det9_time_per_frame     = textdata{57};
        
        pardata.det10_fname             = textdata{58};
        pardata.det10_fnum              = textdata{59};
        pardata.det10_frames_per_file   = textdata{60};
        pardata.det10_time_per_frame    = textdata{61};
        
        %%% LINE 5
        pardata.scaler1_val     = textdata{62};
        pardata.scaler1_units   = textdata{63};
        pardata.scaler2_val     = textdata{64};
        pardata.scaler2_units   = textdata{65};
        pardata.scaler3_val     = textdata{66};
        pardata.scaler3_units   = textdata{67};
        pardata.scaler4_val     = textdata{68};
        pardata.scaler4_units   = textdata{69};
        pardata.scaler5_val     = textdata{70};
        pardata.scaler5_units   = textdata{71};
        pardata.scaler6_val     = textdata{72};
        pardata.scaler6_units   = textdata{73};
        pardata.scaler7_val     = textdata{74};
        pardata.scaler7_units   = textdata{75};
        pardata.scaler8_val     = textdata{76};
        pardata.scaler8_units   = textdata{77};
        pardata.scaler9_val     = textdata{78};
        pardata.scaler9_units   = textdata{79};
        pardata.scaler10_val    = textdata{80};
        pardata.scaler10_units  = textdata{81};
        
        %%% LINE 6
        pardata.samX        = textdata{82};
        pardata.samY        = textdata{83};
        pardata.samZ        = textdata{84};
        pardata.aX          = textdata{85};
        pardata.aY          = textdata{86};
        pardata.aZ          = textdata{87};
        pardata.samX2       = textdata{88};
        pardata.samY2       = textdata{89};
        pardata.samZ2       = textdata{90};
        pardata.samOther    = textdata{91};
        
        %%% LINE 7
        pardata.det1_pos1   = textdata{92};
        pardata.det1_pos2   = textdata{93};
        pardata.det1_pos3   = textdata{94};
        
        pardata.det2_pos1   = textdata{95};
        pardata.det2_pos2   = textdata{96};
        pardata.det2_pos3   = textdata{97};
        
        pardata.det3_pos1   = textdata{98};
        pardata.det3_pos2   = textdata{99};
        pardata.det3_pos3   = textdata{100};
        
        pardata.det4_pos1   = textdata{101};
        pardata.det4_pos2   = textdata{102};
        pardata.det4_pos3   = textdata{103};
        
        pardata.det5_pos1   = textdata{104};
        pardata.det5_pos2   = textdata{105};
        pardata.det5_pos3   = textdata{106};
        
        pardata.det6_pos1   = textdata{107};
        pardata.det6_pos2   = textdata{108};
        pardata.det6_pos3   = textdata{109};
       
        pardata.det7_pos1   = textdata{110};
        pardata.det7_pos2   = textdata{111};
        pardata.det7_pos3   = textdata{112};
        
        pardata.det8_pos1   = textdata{113};
        pardata.det8_pos2   = textdata{114};
        pardata.det8_pos3   = textdata{115};
        
        pardata.det9_pos1   = textdata{116};
        pardata.det9_pos2   = textdata{117};
        pardata.det9_pos3   = textdata{118};
        
        pardata.det10_pos1  = textdata{119};
        pardata.det10_pos2  = textdata{120};
        pardata.det10_pos3  = textdata{121};
        
        %%% LINE 8
        pardata.hex_pos1    = textdata{122};
        pardata.hex_pos2    = textdata{123};
        pardata.hex_pos3    = textdata{124};
        pardata.hex_pos4    = textdata{125};
        pardata.hex_pos5    = textdata{126};
        pardata.hex_pos6    = textdata{127};
        pardata.hex_pos7    = textdata{128};
        
        %%% LINE 9
        pardata.slit1_V_size    = textdata{129};
        pardata.slit1_V_pos     = textdata{130};
        pardata.slit1_H_size    = textdata{131};
        pardata.slit1_H_pos     = textdata{132};
        
        pardata.slit2_V_size    = textdata{133};
        pardata.slit2_V_pos     = textdata{134};
        pardata.slit2_H_size    = textdata{135};
        pardata.slit2_H_pos     = textdata{136};
        
        pardata.slit3_V_size    = textdata{137};
        pardata.slit3_V_pos     = textdata{138};
        pardata.slit3_H_size    = textdata{139};
        pardata.slit3_H_pos     = textdata{140};
        
        pardata.slit4_V_size    = textdata{141};
        pardata.slit4_V_pos     = textdata{142};
        pardata.slit4_H_size    = textdata{143};
        pardata.slit4_H_pos     = textdata{144};
        
        pardata.slit5_V_size    = textdata{145};
        pardata.slit5_V_pos     = textdata{146};
        pardata.slit5_H_size    = textdata{147};
        pardata.slit5_H_pos     = textdata{148};
       
        pardata.slit6_V_size    = textdata{149};
        pardata.slit6_V_pos     = textdata{150};
        pardata.slit6_H_size    = textdata{151};
        pardata.slit6_H_pos     = textdata{152};
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata{153};
        pardata.lens1_pos2  = textdata{154};
        
        pardata.lens2_pos1  = textdata{155};
        pardata.lens2_pos2  = textdata{156};
        
        pardata.lens3_pos1  = textdata{157};
        pardata.lens3_pos2  = textdata{158};
        
        pardata.lens4_pos1  = textdata{159};
        pardata.lens4_pos2  = textdata{160};
        
        %%% LINE 11
        pardata.encoder1    = textdata{161};
        pardata.encoder2    = textdata{162};
        pardata.encoder3    = textdata{163};
        pardata.encoder4    = textdata{164};
        pardata.encoder5    = textdata{165};
        pardata.encoder6    = textdata{166};
        pardata.encoder7    = textdata{167};
        pardata.encoder8    = textdata{168};
        pardata.encoder9    = textdata{169};
        pardata.encoder10   = textdata{170};
        
        %%% LINE 12
        pardata.ev1     = textdata{171};
        pardata.ev2     = textdata{172};
        pardata.ev3     = textdata{173};
        pardata.ev4     = textdata{174};
        pardata.ev5     = textdata{175};
        pardata.ev6     = textdata{176};
        pardata.ev7     = textdata{177};
        pardata.ev8     = textdata{178};
        pardata.ev9     = textdata{179};
        pardata.ev10    = textdata{180};
    case 'mpe_standard'
        fmtstring   = ['%s %s %d %s %d ' ...
            '%u %f %f %f %f %f %f %f ' ....
            '%s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f'];
        
        %%% READ IN DATA USING FORMAT STRING
        fid     = fopen(fname, 'r','n');
        textdata  = textscan(fid, fmtstring);
        fclose(fid);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        pardata.epoch_time  = textdata{6};
        pardata.integ_time  = textdata{7};
        pardata.Iring       = textdata{8};
        pardata.und_gap     = textdata{9};
        pardata.energy      = textdata{10};
        pardata.energy_cal  = textdata{11};
        pardata.foil_pos    = textdata{12};
        pardata.atten_pos   = textdata{13};
        
        pardata.det1_fname              = textdata{14};
        pardata.det1_fnum               = textdata{15};
        pardata.det1_frames_per_file    = textdata{16};
        pardata.det1_time_per_frame     = textdata{17};
        
        pardata.det2_fname              = textdata{18};
        pardata.det2_fnum               = textdata{19};
        pardata.det2_frames_per_file    = textdata{20};
        pardata.det2_time_per_frame     = textdata{21};
        
        pardata.det3_fname              = textdata{22};
        pardata.det3_fnum               = textdata{23};
        pardata.det3_frames_per_file    = textdata{24};
        pardata.det3_time_per_frame     = textdata{25};
        
        pardata.det4_fname              = textdata{26};
        pardata.det4_fnum               = textdata{27};
        pardata.det4_frames_per_file    = textdata{28};
        pardata.det4_time_per_frame     = textdata{29};
        
        pardata.det5_fname              = textdata{30};
        pardata.det5_fnum               = textdata{31};
        pardata.det5_frames_per_file    = textdata{32};
        pardata.det5_time_per_frame     = textdata{33};
        
        pardata.det6_fname              = textdata{34};
        pardata.det6_fnum               = textdata{35};
        pardata.det6_frames_per_file    = textdata{36};
        pardata.det6_time_per_frame     = textdata{37};
        
        pardata.det7_fname              = textdata{38};
        pardata.det7_fnum               = textdata{39};
        pardata.det7_frames_per_file    = textdata{40};
        pardata.det7_time_per_frame     = textdata{41};
        
        pardata.det8_fname              = textdata{42};
        pardata.det8_fnum               = textdata{43};
        pardata.det8_frames_per_file    = textdata{44};
        pardata.det8_time_per_frame     = textdata{45};
        
        pardata.det9_fname              = textdata{46};
        pardata.det9_fnum               = textdata{47};
        pardata.det9_frames_per_file    = textdata{48};
        pardata.det9_time_per_frame     = textdata{49};
        
        pardata.det10_fname             = textdata{50};
        pardata.det10_fnum              = textdata{51};
        pardata.det10_frames_per_file   = textdata{52};
        pardata.det10_time_per_frame    = textdata{53};
        
        pardata.scaler1_val     = textdata{54};
        pardata.scaler1_units   = textdata{55};
        pardata.scaler2_val     = textdata{56};
        pardata.scaler2_units   = textdata{57};
        pardata.scaler3_val     = textdata{58};
        pardata.scaler3_units   = textdata{59};
        pardata.scaler4_val     = textdata{60};
        pardata.scaler4_units   = textdata{61};
        pardata.scaler5_val     = textdata{62};
        pardata.scaler5_units   = textdata{63};
        pardata.scaler6_val     = textdata{64};
        pardata.scaler6_units   = textdata{65};
        pardata.scaler7_val     = textdata{66};
        pardata.scaler7_units   = textdata{67};
        pardata.scaler8_val     = textdata{68};
        pardata.scaler8_units   = textdata{69};
        pardata.scaler9_val     = textdata{70};
        pardata.scaler9_units   = textdata{71};
        pardata.scaler10_val    = textdata{72};
        pardata.scaler10_units  = textdata{73};
        
        pardata.samX        = textdata{74};
        pardata.samY        = textdata{75};
        pardata.samZ        = textdata{76};
        pardata.aX          = textdata{77};
        pardata.aY          = textdata{78};
        pardata.aZ          = textdata{79};
        pardata.samX2       = textdata{80};
        pardata.samY2       = textdata{81};
        pardata.samZ2       = textdata{82};
        pardata.samOther    = textdata{83};
        
        pardata.det1_pos1   = textdata{84};
        pardata.det1_pos2   = textdata{85};
        pardata.det1_pos3   = textdata{86};
        
        pardata.det2_pos1   = textdata{87};
        pardata.det2_pos2   = textdata{88};
        pardata.det2_pos3   = textdata{89};
        
        pardata.det3_pos1   = textdata{90};
        pardata.det3_pos2   = textdata{91};
        pardata.det3_pos3   = textdata{92};
        
        pardata.det4_pos1   = textdata{93};
        pardata.det4_pos2   = textdata{94};
        pardata.det4_pos3   = textdata{95};
        
        pardata.det5_pos1   = textdata{96};
        pardata.det5_pos2   = textdata{97};
        pardata.det5_pos3   = textdata{98};
        
        pardata.det6_pos1   = textdata{99};
        pardata.det6_pos2   = textdata{100};
        pardata.det6_pos3   = textdata{101};
       
        pardata.det7_pos1   = textdata{102};
        pardata.det7_pos2   = textdata{103};
        pardata.det7_pos3   = textdata{104};
        
        pardata.det8_pos1   = textdata{105};
        pardata.det8_pos2   = textdata{106};
        pardata.det8_pos3   = textdata{107};
        
        pardata.det9_pos1   = textdata{108};
        pardata.det9_pos2   = textdata{109};
        pardata.det9_pos3   = textdata{110};
        
        pardata.det10_pos1  = textdata{111};
        pardata.det10_pos2  = textdata{112};
        pardata.det10_pos3  = textdata{113};
        
        pardata.hex_pos1    = textdata{114};
        pardata.hex_pos2    = textdata{115};
        pardata.hex_pos3    = textdata{116};
        pardata.hex_pos4    = textdata{117};
        pardata.hex_pos5    = textdata{118};
        pardata.hex_pos6    = textdata{119};
        pardata.hex_pos7    = textdata{120};
        
        pardata.slit1_V_size    = textdata{121};
        pardata.slit1_V_pos     = textdata{122};
        pardata.slit1_H_size    = textdata{123};
        pardata.slit1_H_pos     = textdata{124};
        
        pardata.slit2_V_size    = textdata{125};
        pardata.slit2_V_pos     = textdata{126};
        pardata.slit2_H_size    = textdata{127};
        pardata.slit2_H_pos     = textdata{128};
        
        pardata.slit3_V_size    = textdata{129};
        pardata.slit3_V_pos     = textdata{130};
        pardata.slit3_H_size    = textdata{131};
        pardata.slit3_H_pos     = textdata{132};
        
        pardata.slit4_V_size    = textdata{133};
        pardata.slit4_V_pos     = textdata{134};
        pardata.slit4_H_size    = textdata{135};
        pardata.slit4_H_pos     = textdata{136};
        
        pardata.slit5_V_size    = textdata{137};
        pardata.slit5_V_pos     = textdata{138};
        pardata.slit5_H_size    = textdata{139};
        pardata.slit5_H_pos     = textdata{140};
        
        pardata.slit6_V_size    = textdata{141};
        pardata.slit6_V_pos     = textdata{142};
        pardata.slit6_H_size    = textdata{143};
        pardata.slit6_H_pos     = textdata{144};
        
        pardata.lens1_pos1  = textdata{145};
        pardata.lens1_pos2  = textdata{146};
        
        pardata.lens2_pos1  = textdata{147};
        pardata.lens2_pos2  = textdata{148};
        
        pardata.lens3_pos1  = textdata{149};
        pardata.lens3_pos2  = textdata{150};
        
        pardata.lens4_pos1  = textdata{151};
        pardata.lens4_pos2  = textdata{152};
        
        pardata.encoder1    = textdata{153};
        pardata.encoder2    = textdata{154};
        pardata.encoder3    = textdata{155};
        pardata.encoder4    = textdata{156};
        pardata.encoder5    = textdata{157};
        pardata.encoder6    = textdata{158};
        pardata.encoder7    = textdata{159};
        pardata.encoder8    = textdata{160};
        pardata.encoder9    = textdata{161};
        pardata.encoder10   = textdata{162};
        
        pardata.ev1     = textdata{163};
        pardata.ev2     = textdata{164};
        pardata.ev3     = textdata{165};
        pardata.ev4     = textdata{166};
        pardata.ev5     = textdata{167};
        pardata.ev6     = textdata{168};
        pardata.ev7     = textdata{169};
        pardata.ev8     = textdata{170};
        pardata.ev9     = textdata{171};
        pardata.ev10    = textdata{172};
    case {'mpe_fastpar_standard', 'mpe_ff_standard'}
        fmtstring   = [ ...
            '%s %s %d %s %d ' ...
            '%s %f %f %f %f %f %f %f %f %f ' ...
            '%s %f %f %d %f ' ...
            '%s %d %d ' ...
            '%f %f %f %f %f ' ...
            '%f %f %f %f %f ' ...
            '%f %f %f %f %f ' ...
            '%f %f %f %f %f ' ...
            '%f %f %f %f %f ' ...
            '%f %f %f %f %f %f ' ...
            '%s %d %f %d %d %d'];
        
        %%% READ IN DATA USING FORMAT STRING
        fid     = fopen(fname, 'r','n');
        textdata  = textscan(fid, fmtstring);
        fclose(fid);
        
        %%% %s %s %8f %8f %8f %8f %8f %8f %8f %8f %8f %6s %5g %5g %4d %5g %12s %05d %05d 
        %12f %12f %12f %12f %12f 
        %12f %12f %12f %8f %8f 
        %8f %8f %12f %12f %12f 
        %12f %12f %12f %12f %12f
        %12f %12f %12f %12f %12f
        %12f %12f %12f %12f %12f
        %8f 
        %s %06d %15.8f %6d %6d %6d
        % enddate, detname, S[SRmA], 
        % A[rze], A[rxe], A[samYE], A[samXE], A[samZE], A[aeroXE], 
        % epics_get("1ide:D1Ch23_calc.VAL"), S[rstr2], motname,\
        % startpos, endpos, OSC["nframes"], OSC["exposure_time"], imgprefix, OSC["first_frame_number"], imgnr,\
        % S[ic1e]/icsec, S[ic2e]/icsec, S[ic3e]/icsec, S[ic4e]/icsec, S[ic5e]/icsec, S[ic6e]/icsec, \
        % S[fedrl], S[fedr2], icsec, displenc, loadcell, stress, tensionmot, S[bposC], S[bposE], \
        % S[bposC], hsizeDS, vsizeDS, hsizeUS, vsizeUS, hposDS, vposDS, hposUS, vposUS, \
        % tiltX, tiltZ, foilwh, attenwh, attenCpos, IF_op, attenXpos, \
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        pardata.Iring       = textdata{7};
        
        pardata.aZ          = textdata{8};
        pardata.aY          = nan;
        pardata.aX          = textdata{9};
        
        pardata.samY        = textdata{10};        
        pardata.samX        = textdata{11};
        pardata.samZ        = textdata{12};
        
        pardata.samX2       = textdata{13};
        pardata.samY2       = nan;
        pardata.samZ2       = nan;
        
        pardata.samOther    = nan;
        
        pardata.scanmtr = textdata{16};
        pardata.scanini = textdata{17};
        pardata.scanfin = textdata{18};
        
        pardata.integ_time  = textdata{20};
        
        pardata.det1_frames_per_file    = textdata{19};
        pardata.det1_time_per_frame     = textdata{20};
        pardata.det1_fname              = textdata{21};
        pardata.det1_fnum               = textdata{22};
        pardata.det2_fname              = nan;
        pardata.det2_fnum               = nan;
        pardata.det2_frames_per_file    = nan;
        pardata.det2_time_per_frame     = textdata{32};
        
        pardata.scaler1_val     = textdata{24};
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = textdata{25};
        pardata.scaler2_units   = nan;
        pardata.scaler3_val     = textdata{26};
        pardata.scaler3_units   = nan;
        pardata.scaler4_val     = textdata{27};
        pardata.scaler4_units   = nan;
        pardata.scaler5_val     = textdata{28};
        pardata.scaler5_units   = nan;
        pardata.scaler6_val     = textdata{29};
        pardata.scaler6_units   = nan;
        pardata.scaler7_val     = nan;
        pardata.scaler7_units   = nan;
        pardata.scaler8_val     = nan;
        pardata.scaler8_units   = nan;
        pardata.scaler9_val     = nan;
        pardata.scaler9_units   = nan;
        pardata.scaler10_val    = nan;
        pardata.scaler10_units  = nan;
        
        pardata.encoder1    = textdata{30};
        pardata.encoder2    = textdata{31};
        pardata.encoder3    = textdata{33};
        pardata.encoder4    = textdata{36};
        pardata.encoder5    = textdata{37};
        pardata.encoder6    = textdata{38};
        pardata.encoder7    = textdata{48};
        pardata.encoder8    = textdata{49};
        pardata.encoder9    = nan;
        pardata.encoder10   = nan;
        
        pardata.ev1     = textdata{14};
        pardata.ev2     = textdata{15};
        pardata.ev3     = textdata{34};
        pardata.ev4     = textdata{35};
        pardata.ev5     = textdata{52};
        pardata.ev6     = textdata{53};
        pardata.ev7     = textdata{54};
        pardata.ev8     = nan;
        pardata.ev9     = nan;
        pardata.ev10    = nan;
        
        pardata.slit1_V_size    = textdata{43};
        pardata.slit1_V_pos     = textdata{47};
        pardata.slit1_H_size    = textdata{42};
        pardata.slit1_H_pos     = textdata{46};
        
        pardata.slit2_V_size    = textdata{41};
        pardata.slit2_V_pos     = textdata{45};
        pardata.slit2_H_size    = textdata{40};
        pardata.slit2_H_pos     = textdata{44};
        
        pardata.slit3_V_size    = nan;
        pardata.slit3_V_pos     = nan;
        pardata.slit3_H_size    = nan;
        pardata.slit3_H_pos     = nan;
        
        pardata.slit4_V_size    = nan;
        pardata.slit4_V_pos     = nan;
        pardata.slit4_H_size    = nan;
        pardata.slit4_H_pos     = nan;
        
        pardata.slit5_V_size    = nan;
        pardata.slit5_V_pos     = nan;
        pardata.slit5_H_size    = nan;
        pardata.slit5_H_pos     = nan;
        
        pardata.slit6_V_size    = nan;
        pardata.slit6_V_pos     = nan;
        pardata.slit6_H_size    = nan;
        pardata.slit6_H_pos     = nan;

        pardata.foil_pos    = textdata{50};
        pardata.atten_pos   = textdata{51};
        
        % sxprefix, sxnum, sxtimestamp, sxcamarraycounter, sximarraycounter, sxfilearraycounter
        pardata.det3_fname              = textdata{55};
        pardata.det3_fnum               = textdata{56};
        pardata.det3_frames_per_file    = nan;
        pardata.det3_time_per_frame     = nan;
        
        pardata.epoch_time  = nan;
        pardata.und_gap     = nan;
        pardata.energy      = nan;
        pardata.energy_cal  = nan;
        
        pardata.det4_fname              = nan;
        pardata.det4_fnum               = nan;
        pardata.det4_frames_per_file    = nan;
        pardata.det4_time_per_frame     = nan;
        
        pardata.det5_fname              = nan;
        pardata.det5_fnum               = nan;
        pardata.det5_frames_per_file    = nan;
        pardata.det5_time_per_frame     = nan;
        
        pardata.det6_fname              = nan;
        pardata.det6_fnum               = nan;
        pardata.det6_frames_per_file    = nan;
        pardata.det6_time_per_frame     = nan;
        
        pardata.det7_fname              = nan;
        pardata.det7_fnum               = nan;
        pardata.det7_frames_per_file    = nan;
        pardata.det7_time_per_frame     = nan;
        
        pardata.det8_fname              = nan;
        pardata.det8_fnum               = nan;
        pardata.det8_frames_per_file    = nan;
        pardata.det8_time_per_frame     = nan;
        
        pardata.det9_fname              = nan;
        pardata.det9_fnum               = nan;
        pardata.det9_frames_per_file    = nan;
        pardata.det9_time_per_frame     = nan;
        
        pardata.det10_fname             = nan;
        pardata.det10_fnum              = nan;
        pardata.det10_frames_per_file   = nan;
        pardata.det10_time_per_frame    = nan;
        
        pardata.det1_pos1   = nan;
        pardata.det1_pos2   = nan;
        pardata.det1_pos3   = nan;
        
        pardata.det2_pos1   = nan;
        pardata.det2_pos2   = nan;
        pardata.det2_pos3   = nan;
        
        pardata.det3_pos1   = nan;
        pardata.det3_pos2   = nan;
        pardata.det3_pos3   = nan;
        
        pardata.det4_pos1   = nan;
        pardata.det4_pos2   = nan;
        pardata.det4_pos3   = nan;
        
        pardata.det5_pos1   = nan;
        pardata.det5_pos2   = nan;
        pardata.det5_pos3   = nan;
        
        pardata.det6_pos1   = nan;
        pardata.det6_pos2   = nan;
        pardata.det6_pos3   = nan;
       
        pardata.det7_pos1   = nan;
        pardata.det7_pos2   = nan;
        pardata.det7_pos3   = nan;
        
        pardata.det8_pos1   = nan;
        pardata.det8_pos2   = nan;
        pardata.det8_pos3   = nan;
        
        pardata.det9_pos1   = nan;
        pardata.det9_pos2   = nan;
        pardata.det9_pos3   = nan;
        
        pardata.det10_pos1  = nan;
        pardata.det10_pos2  = nan;
        pardata.det10_pos3  = nan;
        
        pardata.hex_pos1    = nan;
        pardata.hex_pos2    = nan;
        pardata.hex_pos3    = nan;
        pardata.hex_pos4    = nan;
        pardata.hex_pos5    = nan;
        pardata.hex_pos6    = nan;
        pardata.hex_pos7    = nan;
        
        pardata.lens1_pos1  = nan;
        pardata.lens1_pos2  = nan;
        
        pardata.lens2_pos1  = nan;
        pardata.lens2_pos2  = nan;
        
        pardata.lens3_pos1  = nan;
        pardata.lens3_pos2  = nan;
        
        pardata.lens4_pos1  = nan;
        pardata.lens4_pos2  = nan;
        
        pardata.cal_foil = nan;
    case 'mpe_ff_per_frame'
        fmtstring   = [ ...
            '%s %s %d %s %d ' ...
            '%s %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f ' ...
            '%s %f %f %f %f ' ...
            '%s %d %d %d %f %f ' ...
            '%f %f %f %f'];
        
%         printf("%s %s %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %s %8f %8f %8f %5f %s %05d %05d %04d %12f %12f %12f %12f %12f %15.8f\n",\
%         enddate, 
%         x6detname, x7A[DetX], x8A[DetZ], x9A[NFYE], x10A[bbYE], x11A[ge_x], x12hydraZmot, x13A[imgXE], x14A[imgZE], x15A[imgYE], \
%         x16A[mtsX2E], x17A[rxe], x18A[rze], x19S[fedrl], x20S[fedr2], x21A[mtsYE], x22A[aeroXE], x23A[aero], x24A[samXE], x25A[samZE], \
%         x26motname, x27startpos, x28endpos, x29omegapos, x30OSC["exposure_time"], 
%         x31imgprefix, x32filenum, x33OSC["first_frame_number"], x34iframe+1, x35moncnt[icnt], x36trcnt[icnt], \
%         x37Emoncnt[icnt], x38Etrcnt[icnt], x39cntticks[icnt]/50e6, x40timestamp[icnt])  # We should put in the elapsed time based on the scaler trigger and 10MHz clock (in 1id) or 50 MHz (in 1ide)
        
        %%% READ IN DATA USING FORMAT STRING
        fid     = fopen(fname, 'r','n');
        textdata  = textscan(fid, fmtstring);
        fclose(fid);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        pardata.epoch_time  = textdata{40};
        pardata.integ_time  = textdata{39};
        pardata.Iring       = nan;
        pardata.und_gap     = nan;
        pardata.energy      = nan;
        pardata.energy_cal  = nan;
        pardata.foil_pos    = nan;
        pardata.atten_pos   = nan;
        
        pardata.det1_fname              = textdata{31};
        pardata.det1_fnum               = textdata{32};
        pardata.det1_frames_per_file    = nan;
        pardata.det1_time_per_frame     = textdata{30};
        
        pardata.det2_fname              = nan;
        pardata.det2_fnum               = nan;
        pardata.det2_frames_per_file    = nan;
        pardata.det2_time_per_frame     = nan;
        
        pardata.det3_fname              = nan;
        pardata.det3_fnum               = nan;
        pardata.det3_frames_per_file    = nan;
        pardata.det3_time_per_frame     = nan;
        
        pardata.det4_fname              = nan;
        pardata.det4_fnum               = nan;
        pardata.det4_frames_per_file    = nan;
        pardata.det4_time_per_frame     = nan;
        
        pardata.det5_fname              = nan;
        pardata.det5_fnum               = nan;
        pardata.det5_frames_per_file    = nan;
        pardata.det5_time_per_frame     = nan;
        
        pardata.det6_fname              = nan;
        pardata.det6_fnum               = nan;
        pardata.det6_frames_per_file    = nan;
        pardata.det6_time_per_frame     = nan;
        
        pardata.det7_fname              = nan;
        pardata.det7_fnum               = nan;
        pardata.det7_frames_per_file    = nan;
        pardata.det7_time_per_frame     = nan;
        
        pardata.det8_fname              = nan;
        pardata.det8_fnum               = nan;
        pardata.det8_frames_per_file    = nan;
        pardata.det8_time_per_frame     = nan;
        
        pardata.det9_fname              = nan;
        pardata.det9_fnum               = nan;
        pardata.det9_frames_per_file    = nan;
        pardata.det9_time_per_frame     = nan;
        
        pardata.det10_fname             = nan;
        pardata.det10_fnum              = nan;
        pardata.det10_frames_per_file   = nan;
        pardata.det10_time_per_frame    = nan;
        
        pardata.scaler1_val     = textdata{35};
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = textdata{36};
        pardata.scaler2_units   = nan;
        pardata.scaler3_val     = textdata{37};
        pardata.scaler3_units   = nan;
        pardata.scaler4_val     = textdata{38};
        pardata.scaler4_units   = nan;
        pardata.scaler5_val     = nan;
        pardata.scaler5_units   = nan;
        pardata.scaler6_val     = nan;
        pardata.scaler6_units   = nan;
        pardata.scaler7_val     = nan;
        pardata.scaler7_units   = nan;
        pardata.scaler8_val     = nan;
        pardata.scaler8_units   = nan;
        pardata.scaler9_val     = nan;
        pardata.scaler9_units   = nan;
        pardata.scaler10_val    = nan;
        pardata.scaler10_units  = nan;
        
        pardata.samX        = textdata{24};
        pardata.samY        = textdata{21};
        pardata.samZ        = textdata{25};
        
        pardata.aX          = textdata{17};
        pardata.aY          = textdata{29};
        pardata.aZ          = textdata{18};
        
        pardata.samX2       = textdata{16};
        pardata.samY2       = nan;
        pardata.samZ2       = textdata{23};
        pardata.samOther    = textdata{22};
        
        pardata.scanmtr = textdata{26};
        pardata.scanini = textdata{27};
        pardata.scanfin = textdata{28};
        
        pardata.det1_pos1   = textdata{11};
        pardata.det1_pos2   = nan;
        pardata.det1_pos3   = textdata{12};
        
        pardata.det2_pos1   = textdata{11};
        pardata.det2_pos2   = nan;
        pardata.det2_pos3   = textdata{12};
        
        pardata.det3_pos1   = textdata{11};
        pardata.det3_pos2   = nan;
        pardata.det3_pos3   = textdata{12};
        
        pardata.det4_pos1   = textdata{11};
        pardata.det4_pos2   = nan;
        pardata.det4_pos3   = textdata{12};
        
        pardata.det5_pos1   = textdata{13};
        pardata.det5_pos2   = textdata{15};
        pardata.det5_pos3   = textdata{14};
        
        pardata.det6_pos1   = textdata{7};
        pardata.det6_pos2   = textdata{9};
        pardata.det6_pos3   = textdata{8};
       
        pardata.det7_pos1   = nan;
        pardata.det7_pos2   = nan;
        pardata.det7_pos3   = nan;
        
        pardata.det8_pos1   = nan;
        pardata.det8_pos2   = nan;
        pardata.det8_pos3   = nan;
        
        pardata.det9_pos1   = nan;
        pardata.det9_pos2   = nan;
        pardata.det9_pos3   = nan;
        
        pardata.det10_pos1  = nan;
        pardata.det10_pos2  = nan;
        pardata.det10_pos3  = nan;
        
        pardata.hex_pos1    = nan;
        pardata.hex_pos2    = nan;
        pardata.hex_pos3    = nan;
        pardata.hex_pos4    = nan;
        pardata.hex_pos5    = nan;
        pardata.hex_pos6    = nan;
        pardata.hex_pos7    = nan;
        
        pardata.slit1_V_size    = nan;
        pardata.slit1_V_pos     = nan;
        pardata.slit1_H_size    = nan;
        pardata.slit1_H_pos     = nan;
        
        pardata.slit2_V_size    = nan;
        pardata.slit2_V_pos     = nan;
        pardata.slit2_H_size    = nan;
        pardata.slit2_H_pos     = nan;
        
        pardata.slit3_V_size    = nan;
        pardata.slit3_V_pos     = nan;
        pardata.slit3_H_size    = nan;
        pardata.slit3_H_pos     = nan;
        
        pardata.slit4_V_size    = nan;
        pardata.slit4_V_pos     = nan;
        pardata.slit4_H_size    = nan;
        pardata.slit4_H_pos     = nan;
        
        pardata.slit5_V_size    = nan;
        pardata.slit5_V_pos     = nan;
        pardata.slit5_H_size    = nan;
        pardata.slit5_H_pos     = nan;
        
        pardata.slit6_V_size    = nan;
        pardata.slit6_V_pos     = nan;
        pardata.slit6_H_size    = nan;
        pardata.slit6_H_pos     = nan;
        
        pardata.lens1_pos1  = nan;
        pardata.lens1_pos2  = nan;
        
        pardata.lens2_pos1  = nan;
        pardata.lens2_pos2  = nan;
        
        pardata.lens3_pos1  = nan;
        pardata.lens3_pos2  = nan;
        
        pardata.lens4_pos1  = nan;
        pardata.lens4_pos2  = nan;
        
        pardata.encoder1    = textdata{19};
        pardata.encoder2    = textdata{20};
        pardata.encoder3    = nan;
        pardata.encoder4    = nan;
        pardata.encoder5    = nan;
        pardata.encoder6    = nan;
        pardata.encoder7    = nan;
        pardata.encoder8    = nan;
        pardata.encoder9    = nan;
        pardata.encoder10   = nan;
        
        pardata.ev1     = textdata{6};
        pardata.ev2     = textdata{10};
        pardata.ev3     = textdata{33};
        pardata.ev4     = textdata{24};
        pardata.ev5     = nan;
        pardata.ev6     = nan;
        pardata.ev7     = nan;
        pardata.ev8     = nan;
        pardata.ev9     = nan;
        pardata.ev10    = nan;
        
        pardata.cal_foil = nan;
    
    case 'finegan_dec18'
        %%% READ IN DATA USING FORMAT STRING
        %         textdata    = readtable(fname, 'FileType', 'text', 'Delimiter', ' ');
        %         keyboard
        textdata    = readtable(fname, 'FileType', 'spreadsheet');
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.det_type        = textdata.Var6;
        pardata.det1_fname      = textdata.Var7;
        pardata.det1_fnum      	= textdata.Var8;
        pardata.det1_frames_per_file    = nan;
        pardata.det1_time_per_frame     = textdata.Var18;
        
        pardata.scanmtr = textdata.Var9;
        pardata.scanini = textdata.Var10;
        pardata.scanfin = textdata.Var11;
        
        pardata.samX	= textdata.Var12;
        pardata.samY	= textdata.Var13;
        pardata.samZ	= textdata.Var14;
        
        pardata.aX  = textdata.Var15;
        pardata.aY	= textdata.Var16;
        pardata.aZ  = textdata.Var17;
        
        pardata.scaler1_val     = textdata.Var22;
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = textdata.Var23;
        pardata.scaler2_units   = nan;
        pardata.scaler3_val     = textdata.Var28;
        pardata.scaler3_units   = nan;
        pardata.scaler4_val     = textdata.Var29;
        pardata.scaler4_units   = nan;
        
        pardata.epoch_time  = nan;
        pardata.integ_time  = nan;
        pardata.Iring       = nan;
        pardata.und_gap     = nan;
        pardata.energy      = nan;
        pardata.energy_cal  = nan;
        pardata.foil_pos    = nan;
        pardata.atten_pos   = nan;
        
        pardata.det2_fname              = nan;
        pardata.det2_fnum               = nan;
        pardata.det2_frames_per_file    = nan;
        pardata.det2_time_per_frame     = nan;
        
        pardata.det3_fname              = nan;
        pardata.det3_fnum               = nan;
        pardata.det3_frames_per_file    = nan;
        pardata.det3_time_per_frame     = nan;
        
        pardata.det4_fname              = nan;
        pardata.det4_fnum               = nan;
        pardata.det4_frames_per_file    = nan;
        pardata.det4_time_per_frame     = nan;
        
        pardata.det5_fname              = nan;
        pardata.det5_fnum               = nan;
        pardata.det5_frames_per_file    = nan;
        pardata.det5_time_per_frame     = nan;
        
        pardata.det6_fname              = nan;
        pardata.det6_fnum               = nan;
        pardata.det6_frames_per_file    = nan;
        pardata.det6_time_per_frame     = nan;
        
        pardata.det7_fname              = nan;
        pardata.det7_fnum               = nan;
        pardata.det7_frames_per_file    = nan;
        pardata.det7_time_per_frame     = nan;
        
        pardata.det8_fname              = nan;
        pardata.det8_fnum               = nan;
        pardata.det8_frames_per_file    = nan;
        pardata.det8_time_per_frame     = nan;
        
        pardata.det9_fname              = nan;
        pardata.det9_fnum               = nan;
        pardata.det9_frames_per_file    = nan;
        pardata.det9_time_per_frame     = nan;
        
        pardata.det10_fname             = nan;
        pardata.det10_fnum              = nan;
        pardata.det10_frames_per_file   = nan;
        pardata.det10_time_per_frame    = nan;
        
        pardata.scaler5_val     = nan;
        pardata.scaler5_units   = nan;
        pardata.scaler6_val     = nan;
        pardata.scaler6_units   = nan;
        pardata.scaler7_val     = nan;
        pardata.scaler7_units   = nan;
        pardata.scaler8_val     = nan;
        pardata.scaler8_units   = nan;
        pardata.scaler9_val     = nan;
        pardata.scaler9_units   = nan;
        pardata.scaler10_val    = nan;
        pardata.scaler10_units  = nan;
        
        pardata.samX2       = nan;
        pardata.samY2       = nan;
        pardata.samZ2       = nan;
        pardata.samOther    = nan;
        
        pardata.scanmtr = nan;
        pardata.scanini = nan;
        pardata.scanfin = nan;
        
        pardata.det1_pos1   = nan;
        pardata.det1_pos2   = nan;
        pardata.det1_pos3   = nan;
        
        pardata.det2_pos1   = nan;
        pardata.det2_pos2   = nan;
        pardata.det2_pos3   = nan;
        
        pardata.det3_pos1   = nan;
        pardata.det3_pos2   = nan;
        pardata.det3_pos3   = nan;
        
        pardata.det4_pos1   = nan;
        pardata.det4_pos2   = nan;
        pardata.det4_pos3   = nan;
        
        pardata.det5_pos1   = nan;
        pardata.det5_pos2   = nan;
        pardata.det5_pos3   = nan;
        
        pardata.det6_pos1   = nan;
        pardata.det6_pos2   = nan;
        pardata.det6_pos3   = nan;
       
        pardata.det7_pos1   = nan;
        pardata.det7_pos2   = nan;
        pardata.det7_pos3   = nan;
        
        pardata.det8_pos1   = nan;
        pardata.det8_pos2   = nan;
        pardata.det8_pos3   = nan;
        
        pardata.det9_pos1   = nan;
        pardata.det9_pos2   = nan;
        pardata.det9_pos3   = nan;
        
        pardata.det10_pos1  = nan;
        pardata.det10_pos2  = nan;
        pardata.det10_pos3  = nan;
        
        pardata.hex_pos1    = nan;
        pardata.hex_pos2    = nan;
        pardata.hex_pos3    = nan;
        pardata.hex_pos4    = nan;
        pardata.hex_pos5    = nan;
        pardata.hex_pos6    = nan;
        pardata.hex_pos7    = nan;
        
        pardata.slit1_V_size    = nan;
        pardata.slit1_V_pos     = nan;
        pardata.slit1_H_size    = nan;
        pardata.slit1_H_pos     = nan;
        
        pardata.slit2_V_size    = nan;
        pardata.slit2_V_pos     = nan;
        pardata.slit2_H_size    = nan;
        pardata.slit2_H_pos     = nan;
        
        pardata.slit3_V_size    = nan;
        pardata.slit3_V_pos     = nan;
        pardata.slit3_H_size    = nan;
        pardata.slit3_H_pos     = nan;
        
        pardata.slit4_V_size    = nan;
        pardata.slit4_V_pos     = nan;
        pardata.slit4_H_size    = nan;
        pardata.slit4_H_pos     = nan;
        
        pardata.slit5_V_size    = nan;
        pardata.slit5_V_pos     = nan;
        pardata.slit5_H_size    = nan;
        pardata.slit5_H_pos     = nan;
        
        pardata.slit6_V_size    = nan;
        pardata.slit6_V_pos     = nan;
        pardata.slit6_H_size    = nan;
        pardata.slit6_H_pos     = nan;
        
        pardata.lens1_pos1  = nan;
        pardata.lens1_pos2  = nan;
        
        pardata.lens2_pos1  = nan;
        pardata.lens2_pos2  = nan;
        
        pardata.lens3_pos1  = nan;
        pardata.lens3_pos2  = nan;
        
        pardata.lens4_pos1  = nan;
        pardata.lens4_pos2  = nan;
        
        pardata.encoder1    = nan;
        pardata.encoder2    = nan;
        pardata.encoder3    = nan;
        pardata.encoder4    = nan;
        pardata.encoder5    = nan;
        pardata.encoder6    = nan;
        pardata.encoder7    = nan;
        pardata.encoder8    = nan;
        pardata.encoder9    = nan;
        pardata.encoder10   = nan;
        
        pardata.ev1     = nan;
        pardata.ev2     = nan;
        pardata.ev3     = nan;
        pardata.ev4     = nan;
        pardata.ev5     = nan;
        pardata.ev6     = nan;
        pardata.ev7     = nan;
        pardata.ev8     = nan;
        pardata.ev9     = nan;
        pardata.ev10    = nan;
        
        pardata.cal_foil = nan;
    case 'mli_nov19_c'
        %%% READ IN DATA USING FORMAT STRING
        textdata  = readtable(fname, 'FileType', 'text', 'Delimiter', ' ');
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.det_type        = textdata.Var6;
        pardata.det1_fname      = textdata.Var7;
        pardata.det1_fnum      	= textdata.Var8;
        pardata.det1_frames_per_file    = nan;
        pardata.det1_time_per_frame     = textdata.Var18;
        
        pardata.scanmtr = textdata.Var9;
        pardata.scanini = textdata.Var10;
        pardata.scanfin = textdata.Var11;
        
        pardata.samX	= textdata.Var12;
        pardata.samY	= textdata.Var13;
        pardata.samZ	= textdata.Var14;
        
        pardata.aX  = textdata.Var15;
        pardata.aY	= textdata.Var16;
        pardata.aZ  = textdata.Var17;
        
        pardata.scaler1_val     = textdata.Var22;
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = textdata.Var23;
        pardata.scaler2_units   = nan;
        
        pardata.epoch_time  = nan;
        pardata.integ_time  = nan;
        pardata.Iring       = nan;
        pardata.und_gap     = nan;
        pardata.energy      = nan;
        pardata.energy_cal  = nan;
        pardata.foil_pos    = nan;
        pardata.atten_pos   = nan;
        
        pardata.det2_fname              = nan;
        pardata.det2_fnum               = nan;
        pardata.det2_frames_per_file    = nan;
        pardata.det2_time_per_frame     = nan;
        
        pardata.det3_fname              = nan;
        pardata.det3_fnum               = nan;
        pardata.det3_frames_per_file    = nan;
        pardata.det3_time_per_frame     = nan;
        
        pardata.det4_fname              = nan;
        pardata.det4_fnum               = nan;
        pardata.det4_frames_per_file    = nan;
        pardata.det4_time_per_frame     = nan;
        
        pardata.det5_fname              = nan;
        pardata.det5_fnum               = nan;
        pardata.det5_frames_per_file    = nan;
        pardata.det5_time_per_frame     = nan;
        
        pardata.det6_fname              = nan;
        pardata.det6_fnum               = nan;
        pardata.det6_frames_per_file    = nan;
        pardata.det6_time_per_frame     = nan;
        
        pardata.det7_fname              = nan;
        pardata.det7_fnum               = nan;
        pardata.det7_frames_per_file    = nan;
        pardata.det7_time_per_frame     = nan;
        
        pardata.det8_fname              = nan;
        pardata.det8_fnum               = nan;
        pardata.det8_frames_per_file    = nan;
        pardata.det8_time_per_frame     = nan;
        
        pardata.det9_fname              = nan;
        pardata.det9_fnum               = nan;
        pardata.det9_frames_per_file    = nan;
        pardata.det9_time_per_frame     = nan;
        
        pardata.det10_fname             = nan;
        pardata.det10_fnum              = nan;
        pardata.det10_frames_per_file   = nan;
        pardata.det10_time_per_frame    = nan;
        
        pardata.scaler1_val     = nan;
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = nan;
        pardata.scaler2_units   = nan;
        pardata.scaler3_val     = nan;
        pardata.scaler3_units   = nan;
        pardata.scaler4_val     = nan;
        pardata.scaler4_units   = nan;
        pardata.scaler5_val     = nan;
        pardata.scaler5_units   = nan;
        pardata.scaler6_val     = nan;
        pardata.scaler6_units   = nan;
        pardata.scaler7_val     = nan;
        pardata.scaler7_units   = nan;
        pardata.scaler8_val     = nan;
        pardata.scaler8_units   = nan;
        pardata.scaler9_val     = nan;
        pardata.scaler9_units   = nan;
        pardata.scaler10_val    = nan;
        pardata.scaler10_units  = nan;
        
        pardata.samX2       = nan;
        pardata.samY2       = nan;
        pardata.samZ2       = nan;
        pardata.samOther    = nan;
        
        pardata.scanmtr = nan;
        pardata.scanini = nan;
        pardata.scanfin = nan;
        
        pardata.det1_pos1   = nan;
        pardata.det1_pos2   = nan;
        pardata.det1_pos3   = nan;
        
        pardata.det2_pos1   = nan;
        pardata.det2_pos2   = nan;
        pardata.det2_pos3   = nan;
        
        pardata.det3_pos1   = nan;
        pardata.det3_pos2   = nan;
        pardata.det3_pos3   = nan;
        
        pardata.det4_pos1   = nan;
        pardata.det4_pos2   = nan;
        pardata.det4_pos3   = nan;
        
        pardata.det5_pos1   = nan;
        pardata.det5_pos2   = nan;
        pardata.det5_pos3   = nan;
        
        pardata.det6_pos1   = nan;
        pardata.det6_pos2   = nan;
        pardata.det6_pos3   = nan;
       
        pardata.det7_pos1   = nan;
        pardata.det7_pos2   = nan;
        pardata.det7_pos3   = nan;
        
        pardata.det8_pos1   = nan;
        pardata.det8_pos2   = nan;
        pardata.det8_pos3   = nan;
        
        pardata.det9_pos1   = nan;
        pardata.det9_pos2   = nan;
        pardata.det9_pos3   = nan;
        
        pardata.det10_pos1  = nan;
        pardata.det10_pos2  = nan;
        pardata.det10_pos3  = nan;
        
        pardata.hex_pos1    = nan;
        pardata.hex_pos2    = nan;
        pardata.hex_pos3    = nan;
        pardata.hex_pos4    = nan;
        pardata.hex_pos5    = nan;
        pardata.hex_pos6    = nan;
        pardata.hex_pos7    = nan;
        
        pardata.slit1_V_size    = nan;
        pardata.slit1_V_pos     = nan;
        pardata.slit1_H_size    = nan;
        pardata.slit1_H_pos     = nan;
        
        pardata.slit2_V_size    = nan;
        pardata.slit2_V_pos     = nan;
        pardata.slit2_H_size    = nan;
        pardata.slit2_H_pos     = nan;
        
        pardata.slit3_V_size    = nan;
        pardata.slit3_V_pos     = nan;
        pardata.slit3_H_size    = nan;
        pardata.slit3_H_pos     = nan;
        
        pardata.slit4_V_size    = nan;
        pardata.slit4_V_pos     = nan;
        pardata.slit4_H_size    = nan;
        pardata.slit4_H_pos     = nan;
        
        pardata.slit5_V_size    = nan;
        pardata.slit5_V_pos     = nan;
        pardata.slit5_H_size    = nan;
        pardata.slit5_H_pos     = nan;
        
        pardata.slit6_V_size    = nan;
        pardata.slit6_V_pos     = nan;
        pardata.slit6_H_size    = nan;
        pardata.slit6_H_pos     = nan;
        
        pardata.lens1_pos1  = nan;
        pardata.lens1_pos2  = nan;
        
        pardata.lens2_pos1  = nan;
        pardata.lens2_pos2  = nan;
        
        pardata.lens3_pos1  = nan;
        pardata.lens3_pos2  = nan;
        
        pardata.lens4_pos1  = nan;
        pardata.lens4_pos2  = nan;
        
        pardata.encoder1    = nan;
        pardata.encoder2    = nan;
        pardata.encoder3    = nan;
        pardata.encoder4    = nan;
        pardata.encoder5    = nan;
        pardata.encoder6    = nan;
        pardata.encoder7    = nan;
        pardata.encoder8    = nan;
        pardata.encoder9    = nan;
        pardata.encoder10   = nan;
        
        pardata.ev1     = nan;
        pardata.ev2     = nan;
        pardata.ev3     = nan;
        pardata.ev4     = nan;
        pardata.ev5     = nan;
        pardata.ev6     = nan;
        pardata.ev7     = nan;
        pardata.ev8     = nan;
        pardata.ev9     = nan;
        pardata.ev10    = nan;
        
        pardata.cal_foil = nan;
    otherwise
        disp('format not implemented')
end
