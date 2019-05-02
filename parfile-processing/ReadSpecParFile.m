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
        textdata  = textscan(fid, fmtstring);
        
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
        textdata  = textscan(fid, fmtstring);
        
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
    case 'xu_nov16'
%%s %s %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %6s %5g %5g %4d %5g %12s %05d %05d %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %8f %8f %8f %8f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %4f %4f %8f %8f %s %06d %15.8f %6d %6d %6d        
%       enddate, detname6, S[SRmA]7, A[DetX], A[DetZ], A[NFYE], A[imgXE], A[imgZE], A[imgYE], A[ge_x], hydraZ, A[mtsYE]16, A[mtsX2E]17, A[mtsXE]18, A[aeroXE], A[rxe], A[rze], motname,\
%       startpos23, endpos, OSC["nframes"]25, OSC["exposure_time"]26, imgprefix, OSC["first_frame_number"], imgnr,\
%       S[ic1e]/icsec, S[ic2e]/icsec, S[ic3e]/icsec, S[ic5e]/icsec, S[ic6e]/icsec, S[ic8e]/icsec, S[ic5b]/icsec, S[ic8b]/icsec37,\
%       S[fedrl], S[fedr2]39, icsec, displenc, loadcell, stress, tensionmot, S[bposC]45, S[bposE], \
%       S[bposC], hsizeDS48, vsizeDS, hsizeUS, vsizeUS, hposDS, vposDS, hposUS54, vposUS55, \
%       tiltX, tiltZ, foilwh58, attenwh59, attenCpos, A[bbYE], \
%       sxprefix62, sxnum, sxtimestamp64, sxcamarraycounter, sximarraycounter, sxfilearraycounter
        fmtstring   = ['%s %s %d %s %d %s ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f ' ...
            '%s %f %f ' ...
            '%f %f ' ...
            '%s %f %f ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f ' ...
            '%s %f %u %f %f %f'];
        
        %%% READ IN DATA USING FORMAT STRING
        textdata  = textscan(fid, fmtstring);
        
        dummy   = size(textdata{1},1);
        dummy   = nan(dummy,1);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        pardata.epoch_time  = textdata{64};
        pardata.integ_time  = textdata{26};
        pardata.Iring       = textdata{7};
        pardata.und_gap     = dummy;
        pardata.energy      = dummy;
        pardata.energy_cal  = dummy;
        pardata.foil_pos    = textdata{58};
        pardata.atten_pos   = textdata{59};
        
        pardata.det1_fname              = textdata{27};
        pardata.det1_fnum               = textdata{28};
        pardata.det1_frames_per_file    = textdata{29};
        pardata.det1_time_per_frame     = textdata{26};
        
        pardata.det2_fname              = textdata{27};
        pardata.det2_fnum               = textdata{28};
        pardata.det2_frames_per_file    = textdata{29};
        pardata.det2_time_per_frame     = textdata{26};
        
        pardata.det3_fname              = textdata{27};
        pardata.det3_fnum               = textdata{28};
        pardata.det3_frames_per_file    = textdata{29};
        pardata.det3_time_per_frame     = textdata{26};
        
        pardata.det4_fname              = textdata{27};
        pardata.det4_fnum               = textdata{28};
        pardata.det4_frames_per_file    = textdata{29};
        pardata.det4_time_per_frame     = textdata{26};
        
        pardata.det5_fname              = textdata{62};
        pardata.det5_fnum               = textdata{63};
        pardata.det5_frames_per_file    = 1;
        pardata.det5_time_per_frame     = textdata{26};
        
        pardata.det6_fname              = dummy;
        pardata.det6_fnum               = dummy;
        pardata.det6_frames_per_file    = dummy;
        pardata.det6_time_per_frame     = dummy;
        
        pardata.det7_fname              = dummy;
        pardata.det7_fnum               = dummy;
        pardata.det7_frames_per_file    = dummy;
        pardata.det7_time_per_frame     = dummy;
        
        pardata.det8_fname              = dummy;
        pardata.det8_fnum               = dummy;
        pardata.det8_frames_per_file    = dummy;
        pardata.det8_time_per_frame     = dummy;
        
        pardata.det9_fname              = dummy;
        pardata.det9_fnum               = dummy;
        pardata.det9_frames_per_file    = dummy;
        pardata.det9_time_per_frame     = dummy;
        
        pardata.det10_fname             = dummy;
        pardata.det10_fnum              = dummy;
        pardata.det10_frames_per_file   = dummy;
        pardata.det10_time_per_frame    = dummy;
        
        pardata.scaler1_val     = textdata{30};
        pardata.scaler1_units   = dummy;
        pardata.scaler2_val     = textdata{31};
        pardata.scaler2_units   = dummy;
        pardata.scaler3_val     = textdata{32};
        pardata.scaler3_units   = dummy;
        pardata.scaler4_val     = dummy;
        pardata.scaler4_units   = dummy;
        pardata.scaler5_val     = textdata{33};
        pardata.scaler5_units   = dummy;
        pardata.scaler6_val     = textdata{34};
        pardata.scaler6_units   = dummy;
        pardata.scaler7_val     = dummy;
        pardata.scaler7_units   = dummy;
        pardata.scaler8_val     = textdata{35};
        pardata.scaler8_units   = dummy;
        pardata.scaler9_val     = textdata{36};
        pardata.scaler9_units   = dummy;
        pardata.scaler10_val    = textdata{37};
        pardata.scaler10_units  = dummy;
        
        pardata.samX        = textdata{16};
        pardata.samY        = textdata{17};
        pardata.samZ        = textdata{18};
        pardata.aX          = textdata{19};
        pardata.aY          = textdata{20};
        pardata.aZ          = textdata{21};
        pardata.samX2       = textdata{23};
        pardata.samY2       = dummy;
        pardata.samZ2       = dummy;
        pardata.samOther    = textdata{24};
        
        pardata.det1_pos1   = textdata{14};
        pardata.det1_pos2   = dummy;
        pardata.det1_pos3   = textdata{15};
        
        pardata.det2_pos1   = textdata{14};
        pardata.det2_pos2   = dummy;
        pardata.det2_pos3   = textdata{15};
        
        pardata.det3_pos1   = textdata{14};
        pardata.det3_pos2   = dummy;
        pardata.det3_pos3   = textdata{15};
        
        pardata.det4_pos1   = textdata{14};
        pardata.det4_pos2   = dummy;
        pardata.det4_pos3   = textdata{15};
        
        pardata.det5_pos1   = dummy;
        pardata.det5_pos2   = dummy;
        pardata.det5_pos3   = dummy;
        
        pardata.det6_pos1   = dummy;
        pardata.det6_pos2   = dummy;
        pardata.det6_pos3   = dummy;
       
        pardata.det7_pos1   = dummy;
        pardata.det7_pos2   = dummy;
        pardata.det7_pos3   = dummy;
        
        pardata.det8_pos1   = dummy;
        pardata.det8_pos2   = dummy;
        pardata.det8_pos3   = dummy;
        
        pardata.det9_pos1   = dummy;
        pardata.det9_pos2   = dummy;
        pardata.det9_pos3   = dummy;
        
        pardata.det10_pos1  = dummy;
        pardata.det10_pos2  = dummy;
        pardata.det10_pos3  = dummy;
        
        pardata.hex_pos1    = dummy;
        pardata.hex_pos2    = dummy;
        pardata.hex_pos3    = dummy;
        pardata.hex_pos4    = dummy;
        pardata.hex_pos5    = dummy;
        pardata.hex_pos6    = dummy;
        pardata.hex_pos7    = dummy;
        
        pardata.slit1_V_size    = textdata{51};
        pardata.slit1_V_pos     = textdata{55};
        pardata.slit1_H_size    = textdata{50};
        pardata.slit1_H_pos     = textdata{54};
        
        pardata.slit2_V_size    = textdata{49};
        pardata.slit2_V_pos     = textdata{53};
        pardata.slit2_H_size    = textdata{48};
        pardata.slit2_H_pos     = textdata{52};
        
        pardata.slit3_V_size    = dummy;
        pardata.slit3_V_pos     = dummy;
        pardata.slit3_H_size    = dummy;
        pardata.slit3_H_pos     = dummy;
        
        pardata.slit4_V_size    = dummy;
        pardata.slit4_V_pos     = dummy;
        pardata.slit4_H_size    = dummy;
        pardata.slit4_H_pos     = dummy;
        
        pardata.slit5_V_size    = dummy;
        pardata.slit5_V_pos     = dummy;
        pardata.slit5_H_size    = dummy;
        pardata.slit5_H_pos     = dummy;
        
        pardata.slit6_V_size    = dummy;
        pardata.slit6_V_pos     = dummy;
        pardata.slit6_H_size    = dummy;
        pardata.slit6_H_pos     = dummy;
        
        pardata.lens1_pos1  = dummy;
        pardata.lens1_pos2  = dummy;
        
        pardata.lens2_pos1  = dummy;
        pardata.lens2_pos2  = dummy;
        
        pardata.lens3_pos1  = dummy;
        pardata.lens3_pos2  = dummy;
        
        pardata.lens4_pos1  = dummy;
        pardata.lens4_pos2  = dummy;
        
        pardata.encoder1    = textdata{57};
        pardata.encoder2    = textdata{56};
        pardata.encoder3    = textdata{47};
        pardata.encoder4    = textdata{46};
        pardata.encoder5    = textdata{45};
        pardata.encoder6    = textdata{39};
        pardata.encoder7    = textdata{38};
        pardata.encoder8    = dummy;
        pardata.encoder9    = dummy;
        pardata.encoder10   = dummy;
        
        pardata.ev1     = textdata{41};
        pardata.ev2     = textdata{42};
        pardata.ev3     = textdata{43};
        pardata.ev4     = textdata{44};
        pardata.ev5     = dummy;
        pardata.ev6     = dummy;
        pardata.ev7     = dummy;
        pardata.ev8     = dummy;
        pardata.ev9     = dummy;
        pardata.ev10    = dummy;
    case 'kaoumi_oct17_fastpar'
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
        
        % pardata.det7_fname              = textdata{38};
        pardata.det7_fnum               = textdata{38};
        pardata.det7_frames_per_file    = textdata{39};
        pardata.det7_time_per_frame     = textdata{40};
        
        pardata.det8_fname              = textdata{41};
        pardata.det8_fnum               = textdata{42};
        pardata.det8_frames_per_file    = textdata{43};
        pardata.det8_time_per_frame     = textdata{44};
        
        pardata.det9_fname              = textdata{45};
        pardata.det9_fnum               = textdata{46};
        pardata.det9_frames_per_file    = textdata{47};
        pardata.det9_time_per_frame     = textdata{48};
        
        pardata.det10_fname             = textdata{49};
        pardata.det10_fnum              = textdata{50};
        pardata.det10_frames_per_file   = textdata{51};
        pardata.det10_time_per_frame    = textdata{52};
        
        pardata.scaler1_val     = textdata{53};
        pardata.scaler1_units   = textdata{54};
        pardata.scaler2_val     = textdata{55};
        pardata.scaler2_units   = textdata{56};
        pardata.scaler3_val     = textdata{57};
        pardata.scaler3_units   = textdata{58};
        pardata.scaler4_val     = textdata{59};
        pardata.scaler4_units   = textdata{60};
        pardata.scaler5_val     = textdata{61};
        pardata.scaler5_units   = textdata{62};
        pardata.scaler6_val     = textdata{63};
        pardata.scaler6_units   = textdata{64};
        pardata.scaler7_val     = textdata{65};
        pardata.scaler7_units   = textdata{66};
        pardata.scaler8_val     = textdata{67};
        pardata.scaler8_units   = textdata{68};
        pardata.scaler9_val     = textdata{69};
        pardata.scaler9_units   = textdata{70};
        pardata.scaler10_val    = textdata{71};
        pardata.scaler10_units  = textdata{72};
       
        pardata.samX        = textdata{73};
        pardata.samY        = textdata{74};
        pardata.samZ        = textdata{75};
        pardata.aX          = textdata{76};
        pardata.aY          = textdata{77};
        pardata.aZ          = textdata{78};
        pardata.samX2       = textdata{79};
        pardata.samY2       = textdata{80};
        pardata.samZ2       = textdata{81};
        pardata.samOther    = textdata{82};
       
        pardata.det1_pos1   = textdata{83};
        pardata.det1_pos2   = textdata{84};
        pardata.det1_pos3   = textdata{85};
        
        pardata.det2_pos1   = textdata{86};
        pardata.det2_pos2   = textdata{87};
        pardata.det2_pos3   = textdata{88};
        
        pardata.det3_pos1   = textdata{89};
        pardata.det3_pos2   = textdata{90};
        pardata.det3_pos3   = textdata{91};
        
        pardata.det4_pos1   = textdata{92};
        pardata.det4_pos2   = textdata{93};
        pardata.det4_pos3   = textdata{94};
        
        pardata.det5_pos1   = textdata{95};
        pardata.det5_pos2   = textdata{96};
        pardata.det5_pos3   = textdata{97};
        
        pardata.det6_pos1   = textdata{98};
        pardata.det6_pos2   = textdata{99};
        pardata.det6_pos3   = textdata{100};
       
        pardata.det7_pos1   = textdata{101};
        pardata.det7_pos2   = textdata{102};
        pardata.det7_pos3   = textdata{103};
        
        pardata.det8_pos1   = textdata{104};
        pardata.det8_pos2   = textdata{105};
        pardata.det8_pos3   = textdata{106};
        
        pardata.det9_pos1   = textdata{107};
        pardata.det9_pos2   = textdata{108};
        pardata.det9_pos3   = textdata{109};
        
        pardata.det10_pos1  = textdata{110};
        pardata.det10_pos2  = textdata{111};
        pardata.det10_pos3  = textdata{112};
       
        pardata.hex_pos1    = textdata{113};
        pardata.hex_pos2    = textdata{114};
        pardata.hex_pos3    = textdata{115};
        pardata.hex_pos4    = textdata{116};
        pardata.hex_pos5    = textdata{117};
        pardata.hex_pos6    = textdata{118};
        pardata.hex_pos7    = textdata{119};
        
        pardata.slit1_V_size    = textdata{120};
        pardata.slit1_V_pos     = textdata{121};
        pardata.slit1_H_size    = textdata{122};
        pardata.slit1_H_pos     = textdata{123};
        
        pardata.slit2_V_size    = textdata{124};
        pardata.slit2_V_pos     = textdata{125};
        pardata.slit2_H_size    = textdata{126};
        pardata.slit2_H_pos     = textdata{127};
        
        pardata.slit3_V_size    = textdata{128};
        pardata.slit3_V_pos     = textdata{129};
        pardata.slit3_H_size    = textdata{130};
        pardata.slit3_H_pos     = textdata{131};
        
        pardata.slit4_V_size    = textdata{132};
        pardata.slit4_V_pos     = textdata{133};
        pardata.slit4_H_size    = textdata{134};
        pardata.slit4_H_pos     = textdata{135};
        
        pardata.slit5_V_size    = textdata{136};
        pardata.slit5_V_pos     = textdata{137};
        pardata.slit5_H_size    = textdata{138};
        pardata.slit5_H_pos     = textdata{139};
        
        pardata.slit6_V_size    = textdata{140};
        pardata.slit6_V_pos     = textdata{141};
        pardata.slit6_H_size    = textdata{142};
        pardata.slit6_H_pos     = textdata{143};
        
        pardata.lens1_pos1  = textdata{144};
        pardata.lens1_pos2  = textdata{145};
        
        pardata.lens2_pos1  = textdata{146};
        pardata.lens2_pos2  = textdata{147};
        
        pardata.lens3_pos1  = textdata{148};
        pardata.lens3_pos2  = textdata{149};
        
        pardata.lens4_pos1  = textdata{150};
        pardata.lens4_pos2  = textdata{151};
        
        pardata.encoder1    = textdata{152};
        pardata.encoder2    = textdata{153};
        pardata.encoder3    = textdata{154};
        pardata.encoder4    = textdata{155};
        pardata.encoder5    = textdata{156};
        pardata.encoder6    = textdata{157};
        pardata.encoder7    = textdata{158};
        pardata.encoder8    = textdata{159};
        pardata.encoder9    = textdata{160};
        pardata.encoder10   = textdata{161};
        
        pardata.ev1     = textdata{162};
        pardata.ev2     = textdata{163};
        pardata.ev3     = textdata{164};
        pardata.ev4     = textdata{165};
        pardata.ev5     = textdata{166};
        pardata.ev6     = textdata{167};
        pardata.ev7     = textdata{168};
        pardata.ev8     = textdata{169};
        pardata.ev9     = textdata{170};
        pardata.ev10    = textdata{171};
    case 'white_apr19_ff_per_frame'
        fmtstring   = [ ...
            '%s %s %d %s %d ' ...
            '%s %s %f %s ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f ' ...
            ];
        
%         printf("%s %s %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %8f %s %8f %8f %8f %5f %s %05d %05d %04d %12f %12f %12f %12f %12f %15.8f\n",\
%         enddate, 
%         x6detname, x7A[DetX], x8A[DetZ], x9A[NFYE], x10A[bbYE], x11A[ge_x], x12hydraZmot, x13A[imgXE], x14A[imgZE], x15A[imgYE], \
%         x16A[mtsX2E], x17A[rxe], x18A[rze], x19S[fedrl], x20S[fedr2], x21A[mtsYE], x22A[aeroXE], x23A[aero], x24A[samXE], x25A[samZE], \
%         x26motname, x27startpos, x28endpos, x29omegapos, x30OSC["exposure_time"], 
%         x31imgprefix, x32filenum, x33OSC["first_frame_number"], x34iframe+1, x35moncnt[icnt], x36trcnt[icnt], \
%         x37Emoncnt[icnt], x38Etrcnt[icnt], x39cntticks[icnt]/50e6, x40timestamp[icnt])  # We should put in the elapsed time based on the scaler trigger and 10MHz clock (in 1id) or 50 MHz (in 1ide)
        
        %%% READ IN DATA USING FORMAT STRING
        textdata  = textscan(fid, fmtstring)
        
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
        
        pardata.det1_fname              = textdata{7};
        pardata.det1_fnum               = textdata{8};
        pardata.det1_frames_per_file    = nan;
        pardata.det1_time_per_frame     = textdata{33};
        
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
        
        pardata.scaler1_val     = textdata{36};
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = textdata{37};
        pardata.scaler2_units   = nan;
        pardata.scaler3_val     = textdata{38};
        pardata.scaler3_units   = nan;
        pardata.scaler4_val     = textdata{39};
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
        
        pardata.samX        = textdata{28};
        pardata.samY        = textdata{30};
        pardata.samZ        = textdata{29};
        
        pardata.aX          = textdata{21};
        pardata.aY          = textdata{27};
        pardata.aZ          = textdata{22};
        
        pardata.samX2       = textdata{26};
        pardata.samY2       = textdata{25};
        pardata.samZ2       = textdata{23};
        pardata.samOther    = textdata{20};
        
        pardata.scanmtr = textdata{9};
        pardata.scanini = textdata{10};
        pardata.scanfin = textdata{11};
        
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
    case '6bma_hchoo_feb18'
        fmtstring   = ['%s %s %d %s %d ' ...
            '%s %d %f %f %f %f %f %f %f %f %f %f %f'];
        
        %%% READ IN DATA USING FORMAT STRING
        textdata  = textscan(fid, fmtstring);
        
        %%% PARSE DATA
        pardata.day     = textdata{1};
        pardata.month   = textdata{2};
        pardata.date    = textdata{3};
        pardata.time    = textdata{4};
        pardata.year    = textdata{5};
        
        pardata.det1_fname              = textdata{6};
        pardata.det1_fnum               = textdata{7};
        pardata.det1_exp_time           = textdata{8};
        pardata.dead_time1              = textdata{9};
        pardata.dead_time2              = textdata{10};
        
        pardata.xr                      = textdata{11};
        pardata.yr                      = textdata{12};
        pardata.zr                      = textdata{13};
        pardata.ksamx                   = textdata{14};
        pardata.ksamy                   = textdata{15};
        pardata.ksamz                   = textdata{16};
        pardata.preci                   = textdata{17};
        pardata.enc1                    = textdata{18};
    otherwise
        disp('format not implemented')
end
fclose(fid);
