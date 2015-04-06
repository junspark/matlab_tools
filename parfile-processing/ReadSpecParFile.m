function pardata = ReadSpecParFile(fname, varargin)
% ReadSpecParFile - read par file generated from python
%
%   INPUT:
%
%   fname
%       name of the spec generated par file name
%
%   numchs
%       number of channels recorded in the spec generated par file (excluding the
%       date / time channel)
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
%   'Version'               version of spec metadata file (par file).
%                           example: (coratella_feb15). default value is
%                           'none'. 'mpe_standard' is the new MPE metadata
%                           format for 1-BM-B EDD / 6-BM-A EDD / 1-ID data.

% default options
optcell = {...
    'Version', 'none', ...
    'NumChannels', 58, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

fid     = fopen(fname, 'r','n');
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
        textdata  = textscan(fid, fmtstring);
        
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
    case 'mpe_standard'
        fmtstring   = ['%s %s %s %s %d, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, ' ....
            '%s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, %s, %f, %f, %f, ', ...
            '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ' ...
            '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f'];
        
        %%% READ IN DATA USING FORMAT STRING
        textdata  = textscan(fid, fmtstring);
        
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
    case 'coratella_feb15'
        disp(opts.Version)
        numchs      = 58;
        
        fmtstring   = '%s %s %s %s %d %s';
        for i = 1:1:numchs
            fmtstring   = [fmtstring, ' %f'];
        end
        
        %%% READ IN DATA USING FORMAT STRING
        textdata  = textscan(fid, fmtstring);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        pardata.froot   = textdata{6};
        
        pardata.GE1num  = textdata{7};
        pardata.GE2num  = textdata{8};
        pardata.GE3num  = textdata{9};
        pardata.GE4num  = textdata{10};
        pardata.GE5num  = textdata{11};
        
        pardata.ExpTime = textdata{12};
        pardata.NumFrames       = textdata{13};
        pardata.SAXS_num_ini    = textdata{14};
        pardata.SAXS_num_fin    = textdata{15};
        
        pardata.conXH   = textdata{16};
        pardata.conYH   = textdata{17};
        pardata.conZH   = textdata{18};
        pardata.conUH   = textdata{19};
        pardata.conVH   = textdata{20};
        pardata.conWH   = textdata{21};
        
        pardata.S0  = textdata{22};
        pardata.S1  = textdata{23};
        pardata.S2  = textdata{24};
        pardata.S3  = textdata{25};
        pardata.S5  = textdata{26};
        pardata.S4  = textdata{27};
        pardata.S6  = textdata{28};
        pardata.S8  = textdata{29};
        pardata.S20 = textdata{30};
        pardata.S18 = textdata{31};
        
        pardata.p1Hs    = textdata{32};
        pardata.p1Vs    = textdata{33};
        pardata.p2Hs    = textdata{34};
        pardata.p2Vs    = textdata{35};
        
        pardata.Iring   = textdata{36};
        pardata.Energy  = textdata{37};
        pardata.Energy_cal  = textdata{38};
        
        pardata.preamp1 = textdata{39};
        pardata.preamp2 = textdata{40};
        pardata.preamp3 = textdata{41};
        pardata.preamp4 = textdata{42};
        pardata.preamp5 = textdata{43};
        pardata.preamp6 = textdata{44};
        pardata.preamp7 = textdata{45};
        pardata.preamp8 = textdata{46};
        
        pardata.samXE       = textdata{47};
        pardata.samYE       = textdata{48};
        pardata.samZE       = textdata{49};
        pardata.samRXE      = textdata{50};
        pardata.samRZE      = textdata{51};
        pardata.samX2E      = textdata{52};
        pardata.samZ2E      = textdata{53};
        pardata.phi         = textdata{54};
        
        pardata.keyence1    = textdata{55};
        pardata.keyence2    = textdata{56};
        pardata.cross       = textdata{57};
        pardata.load        = textdata{58};
        pardata.mts3        = textdata{59};
        pardata.mts4        = textdata{60};
        
        pardata.temp1       = textdata{61};
        pardata.temp2       = textdata{62};
        pardata.temp3       = textdata{63};
        pardata.fur_output  = textdata{64};
    case 'gao_mar15'
        disp(opts.Version)
        fmtstring   = '%s %s %s %s %d %s %f %f %f %f %f %f %f %f %f %s %f %f %d %f %s %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f';
        
        %%% READ IN DATA USING FORMAT STRING
        textdata  = textscan(fid, fmtstring);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        pardata.detname = textdata{6};
        pardata.Iring   = textdata{7};
        
        pardata.NFYE    = textdata{8};
        pardata.DetZ    = textdata{9};
        pardata.aeroXE  = textdata{10};
        pardata.samYE   = textdata{11};
        pardata.samXE   = textdata{12};
        pardata.samZE   = textdata{13};
        pardata.vff_r   = textdata{14};
        pardata.vff_eta = textdata{15};
        pardata.scan_mtr    = textdata{16};
        pardata.scan_start  = textdata{17};
        pardata.scan_end    = textdata{18};
        pardata.NumFrames   = textdata{19};
        pardata.ExpTime     = textdata{20};
        pardata.froot       = textdata{21};
        
        pardata.frame_number_start  = textdata{22};
        pardata.GE5num              = textdata{23};
        
        pardata.S1  = textdata{24};
        pardata.S2  = textdata{25};
        pardata.S3  = textdata{26};
        pardata.S4  = textdata{27};
        pardata.S5  = textdata{28};
        pardata.S6  = textdata{29};
        
        pardata.fedrl   = textdata{30};
        pardata.fedr2   = textdata{31};
        pardata.S0      = textdata{32};
        
        pardata.displenc    = textdata{33};
        pardata.loadcell    = textdata{34};
        pardata.stress      = textdata{35};
        pardata.tension_mtr = textdata{36};
       
        pardata.bposC   = textdata{37};
        pardata.bposE   = textdata{38};
        pardata.bposC   = textdata{39};
        pardata.hsizeDS = textdata{40};
        pardata.vsizeDS = textdata{41};
        pardata.hsizeUS = textdata{42};
        pardata.vsizeUS = textdata{43};
        pardata.hposDS  = textdata{44};
        pardata.vposDS  = textdata{45};
        pardata.hposUS  = textdata{46};
        pardata.vposUS  = textdata{47};
        pardata.tiltX   = textdata{48};
        pardata.tiltZ   = textdata{49};
        pardata.att1    = textdata{50};
        pardata.att2    = textdata{51};
        pardata.att3    = textdata{52};
        pardata.att4    = textdata{53};
        
        pardata.attenXpos   = textdata{54};
    otherwise
        disp('format not implemented')
end
fclose(fid);