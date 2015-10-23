function pardata = ReadPythonParFile(fname, varargin)
% ReadPythonParFile - read par file generated from python
%
%   INPUT:
%
%   fname
%       name of the python generated par file name
%
%   Version
%       versions of the pypar file
%       'current' is the current state
%       
%
%   OUTPUT:
%
%   pardata
%       strip chart data in struct arrary.

% default options
optcell = {...
    'Version', 'feb15', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

fid     = fopen(fname, 'r','n');
header  = fgetl(fid);
switch lower(opts.Version)
    case 'pantleon_oct15_bpm'
    %         Date, 
    %         Iring, 
    %         Kohzu DS H size (mm), 
    %         Kohzu DS V size (mm), 
    %         Kohzu DS H pos (mm), 
    %         Kohzu DS V pos (mm), 
    %         Kohzu US H size (mm), 
    %         Kohzu US V size (mm), 
    %         Kohzu US H pos (mm), 
    %         Kohzu US V pos (mm), 
    %         SR BPM H pos (mm), 
    %         SR BPM H angle (urad), 
    %         SR BPM V pos (mm), 
    %         SR BPM V angle (urad), 
    %         XBPM H pos (mm), 
    %         XBPM H angle (urad), 
    %         XBPM V pos (mm), 
    %         XBPM V angle (urad), 
    %         IC-B1 (cts), 
    %         IC-B2 (cts), 
    %         IC-B3 (cts), 
    %         IC-B4 (cts), 
    %         IC-B5 (cts), 
    %         IC-B8 (cts), 
    %         IC-B6 (cts), 
    %         IC-B7 (cts),
    %         world angle (urad), 
    %         E_hrm (keV), 
    %         pzt1 (um), 
    %         pzt2 (um), 
    %         bpm_c (mm), 
    %         bpm_e (mm), 
    %         IC-E3 (cts), 
    %         IC-E5 (cts), 
        fmtstring       = '%s %s %d %s %d, %f,  %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f';
        textdata  = textscan(fid, fmtstring);
        
        pardata.hdr                 = header;
        pardata.day                 = textdata{1};
        pardata.month               = textdata{2};
        pardata.date                = textdata{3};
        pardata.time                = textdata{4};
        pardata.year                = textdata{5};
        pardata.Iring               = textdata{6};
        
        pardata.KDHSize             = textdata{7};
        pardata.KDVSize             = textdata{8};
        pardata.KDHPos              = textdata{9};
        pardata.KDVPos              = textdata{10};
        pardata.KUHSize             = textdata{11};
        pardata.KUVSize             = textdata{12};
        pardata.KUHPos              = textdata{13};
        pardata.KUVPos              = textdata{14};
        
        pardata.SRBPMHpos           = textdata{15};
        pardata.SRBPMHang           = textdata{16};
        pardata.SRBPMVpos           = textdata{17};
        pardata.SRBPMVang           = textdata{18};
        pardata.XBPMHpos            = textdata{19};
        pardata.XBPMHang            = textdata{20};
        pardata.XBPMVpos            = textdata{21};
        pardata.XBPMVang            = textdata{22};

        pardata.ICB1        = textdata{23};
        pardata.ICB2        = textdata{24};
        pardata.ICB3        = textdata{25};
        pardata.ICB4        = textdata{26};
        pardata.ICB5        = textdata{27};
        pardata.ICB8        = textdata{28};
        pardata.ICB6        = textdata{29};
        pardata.ICB7        = textdata{30};
        
        pardata.world_angle = textdata{31};
        pardata.E_hrm       = textdata{32};
        pardata.pzt1        = textdata{33};
        pardata.pzt2        = textdata{34};
        pardata.bpm_C       = textdata{35};
        pardata.bpm_E       = textdata{36};
        pardata.ICE3        = textdata{37};
        pardata.ICE5        = textdata{38};
    case 'feb15'
        fmtstring       = '%s %s %d %s %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f';
        textdata  = textscan(fid, fmtstring);
        
        pardata.hdr                 = header;
        pardata.day                 = textdata{1};
        pardata.month               = textdata{2};
        pardata.date                = textdata{3};
        pardata.time                = textdata{4};
        pardata.year                = textdata{5};
        pardata.Iring               = textdata{6};
        pardata.ElapsedTime         = textdata{7};
        pardata.ic0                 = textdata{8};
        pardata.ic1                 = textdata{9};
        pardata.NULL                = textdata{10};
        pardata.HorzInstDeadTime    = textdata{11};
        pardata.HorzAveDeadTime     = textdata{12};
        pardata.HorzRealElapTime    = textdata{13};
        pardata.HorzLiveElapTime    = textdata{14};
        pardata.VertInstDeadTime    = textdata{15};
        pardata.VertAveDeadTime     = textdata{16};
        pardata.VertRealElapTime    = textdata{17};
        pardata.VertLiveElapTime    = textdata{18};
        pardata.IterationNumber     = textdata{19};
        pardata.PositionNumber      = textdata{20};
        pardata.FileIDNumber        = textdata{21};
        pardata.HorzFileName        = textdata{22};
        pardata.VertFileName        = textdata{23};
        pardata.samX                = textdata{24};
        pardata.samY                = textdata{25};
        pardata.samY2               = textdata{26};
        pardata.samZ                = textdata{27};
        pardata.HorzSlitB           = textdata{28};
        pardata.HorzSlitT           = textdata{29};
        pardata.HorzSlitI           = textdata{30};
        pardata.HorzSlitO           = textdata{31};
        pardata.HorzTh              = textdata{32};
        pardata.VertSlitB           = textdata{33};
        pardata.VertSlitT           = textdata{34};
        pardata.VertSlitI           = textdata{35};
        pardata.VertSlitO           = textdata{36};
        pardata.VertTh              = textdata{37};
    case 'dec14'
        fmtstring       = '%s %s %d %s %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s %s %f, %f, %f, %f, %f';
        textdata  = textscan(fid, fmtstring);
        
        pardata.hdr                 = header;
        pardata.day                 = textdata{1};
        pardata.month               = textdata{2};
        pardata.date                = textdata{3};
        pardata.time                = textdata{4};
        pardata.year                = textdata{5};
        pardata.Iring               = textdata{6};
        pardata.ElapsedTime         = textdata{7};
        pardata.ic0                 = textdata{8};
        pardata.ic1                 = textdata{9};
        pardata.NULL                = textdata{10};
        pardata.HorzInstDeadTime    = textdata{11};
        pardata.HorzAveDeadTime     = textdata{12};
        pardata.HorzRealElapTime    = textdata{13};
        pardata.HorzLiveElapTime    = textdata{14};
        pardata.samX                = textdata{15};
        pardata.samY                = textdata{16};
        pardata.samY2               = textdata{17};
        pardata.samZ                = textdata{18};
        pardata.HorzSlitB           = textdata{19};
        pardata.HorzSlitT           = textdata{20};
        pardata.HorzSlitI           = textdata{21};
        pardata.HorzSlitO           = textdata{22};
        pardata.HorzTh              = textdata{23};
        pardata.VertSlitB           = textdata{24};
        pardata.VertSlitT           = textdata{25};
        pardata.VertSlitI           = textdata{26};
        pardata.VertSlitO           = textdata{27};
        pardata.VertTh              = textdata{28};
    case 'mpe_standard'
        fmtstring   = ['%s %s %d %s %d, ' ...
            '%u, %f, %f, %f, %f, %f, %f, %f, ' ....
            '%s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, %s %f, %f, %f, ' ...
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
    otherwise
        disp('format not implemented')
end
fclose(fid);