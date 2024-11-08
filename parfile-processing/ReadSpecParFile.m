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
    
    case 'saxs_waxs_fmt_fastpar_brown_nov19'
        %%% READ IN DATA
        % textdata    = readtable(fname, 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'split');
        opts        = detectImportOptions(fname,'FileType','text');
        textdata    = readtable(fname, opts);
        % pardata     = textdata;
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_fnum               = textdata.Var42;
        pardata.det5_frames_per_file    = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
                    
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
        pardata.det10_fname             = textdata.Var61;
        pardata.det10_fnum              = textdata.Var62;
        pardata.det10_frames_per_file   = textdata.Var63;
        pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var65;
        pardata.scaler1_units   = textdata.Var66;
        pardata.scaler2_val     = textdata.Var67;
        pardata.scaler2_units   = textdata.Var68;
        pardata.scaler3_val     = textdata.Var69;
        pardata.scaler3_units   = textdata.Var70;
        pardata.scaler4_val     = textdata.Var71;
        pardata.scaler4_units   = textdata.Var72;
        pardata.scaler5_val     = textdata.Var73;
        pardata.scaler5_units   = textdata.Var74;
        pardata.scaler6_val     = textdata.Var75;
        pardata.scaler6_units   = textdata.Var76;
        pardata.scaler7_val     = textdata.Var77;
        pardata.scaler7_units   = textdata.Var78;
        pardata.scaler8_val     = textdata.Var79;
        pardata.scaler8_units   = textdata.Var80;
        pardata.scaler9_val     = textdata.Var81;
        pardata.scaler9_units   = textdata.Var82;
        pardata.scaler10_val    = textdata.Var83;
        pardata.scaler10_units  = textdata.Var84;
        pardata.scaler11_val    = textdata.Var85;
        pardata.scaler11_units  = textdata.Var86;
        pardata.scaler12_val    = textdata.Var87;
        pardata.scaler12_units  = textdata.Var88;
        
        pardata.samX        = textdata.Var89;
        pardata.samY        = textdata.Var90;
        pardata.samZ        = textdata.Var91;
        pardata.aX          = textdata.Var92;
        pardata.aY          = textdata.Var93;
        pardata.aZ          = textdata.Var94;
        pardata.samX2       = textdata.Var95;
        pardata.samY2       = textdata.Var96;
        pardata.samZ2       = textdata.Var97;
        pardata.samOther    = textdata.Var98;
        
        pardata.det1_pos1   = textdata.Var99;
        pardata.det1_pos2   = textdata.Var100;
        pardata.det1_pos3   = textdata.Var101;
       
        pardata.det2_pos1   = textdata.Var102;
        pardata.det2_pos2   = textdata.Var103;
        pardata.det2_pos3   = textdata.Var104;
        
        pardata.det3_pos1   = textdata.Var105;
        pardata.det3_pos2   = textdata.Var106;
        pardata.det3_pos3   = textdata.Var107;
        
        pardata.det4_pos1   = textdata.Var108;
        pardata.det4_pos2   = textdata.Var109;
        pardata.det4_pos3   = textdata.Var110;
        
        pardata.det5_pos1   = textdata.Var111;
        pardata.det5_pos2   = textdata.Var112;
        pardata.det5_pos3   = textdata.Var113;
        
        pardata.det6_pos1   = textdata.Var114;
        pardata.det6_pos2   = textdata.Var115;
        pardata.det6_pos3   = textdata.Var116;
      
        pardata.det7_pos1   = textdata.Var117;
        pardata.det7_pos2   = textdata.Var118;
        pardata.det7_pos3   = textdata.Var119;
       
        pardata.det8_pos1   = textdata.Var120;
        pardata.det8_pos2   = textdata.Var121;
        pardata.det8_pos3   = textdata.Var122;
       
        pardata.det9_pos1   = textdata.Var123;
        pardata.det9_pos2   = textdata.Var124;
        pardata.det9_pos3   = textdata.Var125;
        
        pardata.det10_pos1  = textdata.Var126;
        pardata.det10_pos2  = textdata.Var127;
        pardata.det10_pos3  = textdata.Var128;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var129;
        pardata.hex_pos2    = textdata.Var130;
        pardata.hex_pos3    = textdata.Var131;
        pardata.hex_pos4    = textdata.Var132;
        pardata.hex_pos5    = textdata.Var133;
        pardata.hex_pos6    = textdata.Var134;
        pardata.hex_pos7    = textdata.Var135;
        
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var136;
        pardata.slit1_V_pos     = textdata.Var137;
        pardata.slit1_H_size    = textdata.Var138;
        pardata.slit1_H_pos     = textdata.Var139;
       
        pardata.slit2_V_size    = textdata.Var140;
        pardata.slit2_V_pos     = textdata.Var141;
        pardata.slit2_H_size    = textdata.Var142;
        pardata.slit2_H_pos     = textdata.Var143;
        
        pardata.slit3_V_size    = textdata.Var144;
        pardata.slit3_V_pos     = textdata.Var145;
        pardata.slit3_H_size    = textdata.Var146;
        pardata.slit3_H_pos     = textdata.Var147;
        
        pardata.slit4_V_size    = textdata.Var148;
        pardata.slit4_V_pos     = textdata.Var149;
        pardata.slit4_H_size    = textdata.Var150;
        pardata.slit4_H_pos     = textdata.Var151;
        
        pardata.slit5_V_size    = textdata.Var152;
        pardata.slit5_V_pos     = textdata.Var153;
        pardata.slit5_H_size    = textdata.Var154;
        pardata.slit5_H_pos     = textdata.Var155;
        
        pardata.slit6_V_size    = textdata.Var156;
        pardata.slit6_V_pos     = textdata.Var157;
        pardata.slit6_H_size    = textdata.Var158;
        pardata.slit6_H_pos     = textdata.Var159;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var160;
        pardata.lens1_pos2  = textdata.Var161;
       
        pardata.lens2_pos1  = textdata.Var162;
        pardata.lens2_pos2  = textdata.Var163;
       
        pardata.lens3_pos1  = textdata.Var164;
        pardata.lens3_pos2  = textdata.Var165;
       
        pardata.lens4_pos1  = textdata.Var166;
        pardata.lens4_pos2  = textdata.Var167;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var168;
        pardata.encoder2    = textdata.Var169;
        pardata.encoder3    = textdata.Var170;
        % pardata.encoder4    = textdata.Var171;
        pardata.encoder5    = textdata.Var171;
        pardata.encoder6    = textdata.Var172;
        pardata.encoder7    = textdata.Var173;
        pardata.encoder8    = textdata.Var174;
        pardata.encoder9    = textdata.Var175;
        pardata.encoder10   = textdata.Var176;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var177;
        pardata.ev2     = textdata.Var178;
        pardata.ev3     = textdata.Var179;
        pardata.ev4     = textdata.Var180;
        pardata.ev5     = textdata.Var181;
        pardata.ev6     = textdata.Var182;
        pardata.ev7     = textdata.Var183;
        pardata.ev8     = textdata.Var184;
        pardata.ev9     = textdata.Var185;
        pardata.ev10    = textdata.Var186;
    case 'saxs_waxs_fmt_fastpar_v2'
        %%% STARTING park_aug20 (maybe even before when slave motors for fs was enabled)
        
        %%% READ IN DATA
        % textdata    = readtable(fname, 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'split');
        opts        = detectImportOptions(fname,'FileType','text');
        textdata    = readtable(fname, opts);
        
        % pardata = textdata;
        % return
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_fnum               = textdata.Var42;
        pardata.det5_frames_per_file    = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
        
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
        pardata.det10_fname             = textdata.Var61;
        pardata.det10_fnum              = textdata.Var62;
        pardata.det10_frames_per_file   = textdata.Var63;
        pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var65;
        pardata.scaler1_units   = textdata.Var66;
        pardata.scaler2_val     = textdata.Var67;
        pardata.scaler2_units   = textdata.Var68;
        pardata.scaler3_val     = textdata.Var69;
        pardata.scaler3_units   = textdata.Var70;
        pardata.scaler4_val     = textdata.Var71;
        pardata.scaler4_units   = textdata.Var72;
        pardata.scaler5_val     = textdata.Var73;
        pardata.scaler5_units   = textdata.Var74;
        pardata.scaler6_val     = textdata.Var75;
        pardata.scaler6_units   = textdata.Var76;
        pardata.scaler7_val     = textdata.Var77;
        pardata.scaler7_units   = textdata.Var78;
        pardata.scaler8_val     = textdata.Var79;
        pardata.scaler8_units   = textdata.Var80;
        pardata.scaler9_val     = textdata.Var81;
        pardata.scaler9_units   = textdata.Var82;
        pardata.scaler10_val    = textdata.Var83;
        pardata.scaler10_units  = textdata.Var84;
        
        pardata.samX        = textdata.Var85;
        pardata.samY        = textdata.Var86;
        pardata.samZ        = textdata.Var87;
        pardata.aX          = textdata.Var88;
        pardata.aY          = textdata.Var89;
        pardata.aZ          = textdata.Var90;
        pardata.samX2       = textdata.Var91;
        pardata.samY2       = textdata.Var92;
        pardata.samZ2       = textdata.Var93;
        pardata.samOther    = textdata.Var94;
        
        pardata.det1_pos1   = textdata.Var95;
        pardata.det1_pos2   = textdata.Var96;
        pardata.det1_pos3   = textdata.Var97;
       
        pardata.det2_pos1   = textdata.Var98;
        pardata.det2_pos2   = textdata.Var99;
        pardata.det2_pos3   = textdata.Var100;
        
        pardata.det3_pos1   = textdata.Var101;
        pardata.det3_pos2   = textdata.Var102;
        pardata.det3_pos3   = textdata.Var103;
        
        pardata.det4_pos1   = textdata.Var104;
        pardata.det4_pos2   = textdata.Var105;
        pardata.det4_pos3   = textdata.Var106;
        
        pardata.det5_pos1   = textdata.Var107;
        pardata.det5_pos2   = textdata.Var108;
        pardata.det5_pos3   = textdata.Var109;
        
        pardata.det6_pos1   = textdata.Var110;
        pardata.det6_pos2   = textdata.Var111;
        pardata.det6_pos3   = textdata.Var112;
      
        pardata.det7_pos1   = textdata.Var113;
        pardata.det7_pos2   = textdata.Var114;
        pardata.det7_pos3   = textdata.Var115;
       
        pardata.det8_pos1   = textdata.Var116;
        pardata.det8_pos2   = textdata.Var117;
        pardata.det8_pos3   = textdata.Var118;
       
        pardata.det9_pos1   = textdata.Var119;
        pardata.det9_pos2   = textdata.Var120;
        pardata.det9_pos3   = textdata.Var121;
        
        pardata.det10_pos1  = textdata.Var122;
        pardata.det10_pos2  = textdata.Var123;
        pardata.det10_pos3  = textdata.Var124;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var125;
        pardata.hex_pos2    = textdata.Var126;
        pardata.hex_pos3    = textdata.Var127;
        pardata.hex_pos4    = textdata.Var128;
        pardata.hex_pos5    = textdata.Var129;
        pardata.hex_pos6    = textdata.Var130;
        pardata.hex_pos7    = textdata.Var131;
        
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var132;
        pardata.slit1_V_pos     = textdata.Var133;
        pardata.slit1_H_size    = textdata.Var134;
        pardata.slit1_H_pos     = textdata.Var135;
       
        pardata.slit2_V_size    = textdata.Var136;
        pardata.slit2_V_pos     = textdata.Var137;
        pardata.slit2_H_size    = textdata.Var138;
        pardata.slit2_H_pos     = textdata.Var139;
        
        pardata.slit3_V_size    = textdata.Var140;
        pardata.slit3_V_pos     = textdata.Var141;
        pardata.slit3_H_size    = textdata.Var142;
        pardata.slit3_H_pos     = textdata.Var143;
        
        pardata.slit4_V_size    = textdata.Var144;
        pardata.slit4_V_pos     = textdata.Var145;
        pardata.slit4_H_size    = textdata.Var146;
        pardata.slit4_H_pos     = textdata.Var147;
        
        pardata.slit5_V_size    = textdata.Var148;
        pardata.slit5_V_pos     = textdata.Var149;
        pardata.slit5_H_size    = textdata.Var150;
        pardata.slit5_H_pos     = textdata.Var151;
       
        pardata.slit6_V_size    = textdata.Var152;
        pardata.slit6_V_pos     = textdata.Var153;
        pardata.slit6_H_size    = textdata.Var154;
        pardata.slit6_H_pos     = textdata.Var155;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var156;
        pardata.lens1_pos2  = textdata.Var157;
       
        pardata.lens2_pos1  = textdata.Var158;
        pardata.lens2_pos2  = textdata.Var159;
       
        pardata.lens3_pos1  = textdata.Var160;
        pardata.lens3_pos2  = textdata.Var161;
       
        pardata.lens4_pos1  = textdata.Var162;
        pardata.lens4_pos2  = textdata.Var163;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var164;
        pardata.encoder2    = textdata.Var165;
        pardata.encoder3    = textdata.Var166;
        pardata.encoder4    = textdata.Var167;
        pardata.encoder5    = textdata.Var168;
        pardata.encoder6    = textdata.Var169;
        pardata.encoder7    = textdata.Var170;
        pardata.encoder8    = textdata.Var171;
        pardata.encoder9    = textdata.Var172;
        pardata.encoder10   = textdata.Var173;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var174;
        pardata.ev2     = textdata.Var175;
        pardata.ev3     = textdata.Var176;
        pardata.ev4     = textdata.Var177;
        pardata.ev5     = textdata.Var178;
        pardata.ev6     = textdata.Var179;
        pardata.ev7     = textdata.Var180;
        pardata.ev8     = textdata.Var181;
        pardata.ev9     = textdata.Var182;
        pardata.ev10    = textdata.Var183;   
    case 'saxs_waxs_fmt_fastpar_v3'
        %%% STARTING mpe1_oct20
        %%% READ IN DATA
        opts        = detectImportOptions(fname,'FileType','text');
        textdata    = readtable(fname, opts);
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_fnum               = textdata.Var42;
        pardata.det5_frames_per_file    = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
        
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
        pardata.det10_fname             = textdata.Var61;
        pardata.det10_fnum              = textdata.Var62;
        pardata.det10_frames_per_file   = textdata.Var63;
        pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var65;
        pardata.scaler1_units   = textdata.Var66;
        pardata.scaler2_val     = textdata.Var67;
        pardata.scaler2_units   = textdata.Var68;
        pardata.scaler3_val     = textdata.Var69;
        pardata.scaler3_units   = textdata.Var70;
        pardata.scaler4_val     = textdata.Var71;
        pardata.scaler4_units   = textdata.Var72;
        pardata.scaler5_val     = textdata.Var73;
        pardata.scaler5_units   = textdata.Var74;
        pardata.scaler6_val     = textdata.Var75;
        pardata.scaler6_units   = textdata.Var76;
        pardata.scaler7_val     = textdata.Var77;
        pardata.scaler7_units   = textdata.Var78;
        pardata.scaler8_val     = textdata.Var79;
        pardata.scaler8_units   = textdata.Var80;
        pardata.scaler9_val     = textdata.Var81;
        pardata.scaler9_units   = textdata.Var82;
        pardata.scaler10_val    = textdata.Var83;
        pardata.scaler10_units  = textdata.Var84;
        
        pardata.samX        = textdata.Var85;
        pardata.samY        = textdata.Var86;
        pardata.samZ        = textdata.Var87;
        pardata.aX          = textdata.Var88;
        pardata.aY          = textdata.Var89;
        pardata.aZ          = textdata.Var90;
        pardata.samX2       = textdata.Var91;
        pardata.samY2       = textdata.Var92;
        pardata.samZ2       = textdata.Var93;
        pardata.samOther    = textdata.Var94;
        
        pardata.det1_pos1   = textdata.Var95;
        pardata.det1_pos2   = textdata.Var96;
        pardata.det1_pos3   = textdata.Var97;
       
        pardata.det2_pos1   = textdata.Var98;
        pardata.det2_pos2   = textdata.Var99;
        pardata.det2_pos3   = textdata.Var100;
        
        pardata.det3_pos1   = textdata.Var101;
        pardata.det3_pos2   = textdata.Var102;
        pardata.det3_pos3   = textdata.Var103;
        
        pardata.det4_pos1   = textdata.Var104;
        pardata.det4_pos2   = textdata.Var105;
        pardata.det4_pos3   = textdata.Var106;
        
        pardata.det5_pos1   = textdata.Var107;
        pardata.det5_pos2   = textdata.Var108;
        pardata.det5_pos3   = textdata.Var109;
        
        pardata.det6_pos1   = textdata.Var110;
        pardata.det6_pos2   = textdata.Var111;
        pardata.det6_pos3   = textdata.Var112;
      
        pardata.det7_pos1   = textdata.Var113;
        pardata.det7_pos2   = textdata.Var114;
        pardata.det7_pos3   = textdata.Var115;
       
        pardata.det8_pos1   = textdata.Var116;
        pardata.det8_pos2   = textdata.Var117;
        pardata.det8_pos3   = textdata.Var118;
       
        pardata.det9_pos1   = textdata.Var119;
        pardata.det9_pos2   = textdata.Var120;
        pardata.det9_pos3   = textdata.Var121;
        
        pardata.det10_pos1  = textdata.Var122;
        pardata.det10_pos2  = textdata.Var123;
        pardata.det10_pos3  = textdata.Var124;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var125;
        pardata.hex_pos2    = textdata.Var126;
        pardata.hex_pos3    = textdata.Var127;
        pardata.hex_pos4    = textdata.Var128;
        pardata.hex_pos5    = textdata.Var129;
        pardata.hex_pos6    = textdata.Var130;
        pardata.hex_pos7    = textdata.Var131;
        
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var132;
        pardata.slit1_V_pos     = textdata.Var133;
        pardata.slit1_H_size    = textdata.Var134;
        pardata.slit1_H_pos     = textdata.Var135;
       
        pardata.slit2_V_size    = textdata.Var136;
        pardata.slit2_V_pos     = textdata.Var137;
        pardata.slit2_H_size    = textdata.Var138;
        pardata.slit2_H_pos     = textdata.Var139;
        
        pardata.slit3_V_size    = textdata.Var140;
        pardata.slit3_V_pos     = textdata.Var141;
        pardata.slit3_H_size    = textdata.Var142;
        pardata.slit3_H_pos     = textdata.Var143;
        
        pardata.slit4_V_size    = textdata.Var144;
        pardata.slit4_V_pos     = textdata.Var145;
        pardata.slit4_H_size    = textdata.Var146;
        pardata.slit4_H_pos     = textdata.Var147;
        
        pardata.slit5_V_size    = textdata.Var148;
        pardata.slit5_V_pos     = textdata.Var149;
        pardata.slit5_H_size    = textdata.Var150;
        pardata.slit5_H_pos     = textdata.Var151;
       
        pardata.slit6_V_size    = textdata.Var152;
        pardata.slit6_V_pos     = textdata.Var153;
        pardata.slit6_H_size    = textdata.Var154;
        pardata.slit6_H_pos     = textdata.Var155;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var156;
        pardata.lens1_pos2  = textdata.Var157;
       
        pardata.lens2_pos1  = textdata.Var158;
        pardata.lens2_pos2  = textdata.Var159;
       
        pardata.lens3_pos1  = textdata.Var160;
        pardata.lens3_pos2  = textdata.Var161;
       
        pardata.lens4_pos1  = textdata.Var162;
        pardata.lens4_pos2  = textdata.Var163;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var164;
        pardata.encoder2    = textdata.Var165;
        pardata.encoder3    = textdata.Var166;
        pardata.encoder4    = textdata.Var167;
        pardata.encoder5    = textdata.Var168;
        pardata.encoder6    = textdata.Var169;
        pardata.encoder7    = textdata.Var170;
        pardata.encoder8    = textdata.Var171;
        pardata.encoder9    = textdata.Var172;
        pardata.encoder10   = textdata.Var173;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var174;
        pardata.ev2     = textdata.Var175;
        pardata.ev3     = textdata.Var176;
        pardata.ev4     = textdata.Var177;
        pardata.ev5     = textdata.Var178;
        pardata.ev6     = textdata.Var179;
        pardata.ev7     = textdata.Var180;
        pardata.ev8     = textdata.Var181;
        pardata.ev9     = textdata.Var182;
        pardata.ev10    = textdata.Var183;
    case 'saxs_waxs_fmt_fastpar_v4'
        %%% STARTING brown_jun21 (after edits)
        %%% READ IN DATA
        opts        = detectImportOptions(fname,'FileType','text');
        textdata    = readtable(fname, opts);
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_frames_per_file    = textdata.Var42;
        pardata.det5_fnum               = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
        
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
        pardata.det10_fname             = textdata.Var61;
        pardata.det10_fnum              = textdata.Var62;
        pardata.det10_frames_per_file   = textdata.Var63;
        pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var65;
        pardata.scaler1_units   = textdata.Var66;
        pardata.scaler2_val     = textdata.Var67;
        pardata.scaler2_units   = textdata.Var68;
        pardata.scaler3_val     = textdata.Var69;
        pardata.scaler3_units   = textdata.Var70;
        pardata.scaler4_val     = textdata.Var71;
        pardata.scaler4_units   = textdata.Var72;
        pardata.scaler5_val     = textdata.Var73;
        pardata.scaler5_units   = textdata.Var74;
        pardata.scaler6_val     = textdata.Var75;
        pardata.scaler6_units   = textdata.Var76;
        pardata.scaler7_val     = textdata.Var77;
        pardata.scaler7_units   = textdata.Var78;
        pardata.scaler8_val     = textdata.Var79;
        pardata.scaler8_units   = textdata.Var80;
        pardata.scaler9_val     = textdata.Var81;
        pardata.scaler9_units   = textdata.Var82;
        pardata.scaler10_val    = textdata.Var83;
        pardata.scaler10_units  = textdata.Var84;
        pardata.scaler11_val    = textdata.Var85;
        pardata.scaler11_units  = textdata.Var86;
        pardata.scaler12_val    = textdata.Var87;
        pardata.scaler12_units  = textdata.Var88;
        
        pardata.samX        = textdata.Var89;
        pardata.samY        = textdata.Var90;
        pardata.samZ        = textdata.Var91;
        pardata.aX          = textdata.Var92;
        pardata.aY          = textdata.Var93;
        pardata.aZ          = textdata.Var94;
        pardata.samX2       = textdata.Var95;
        pardata.samY2       = textdata.Var96;
        pardata.samZ2       = textdata.Var97;
        pardata.samOther    = textdata.Var98;
        
        pardata.det1_pos1   = textdata.Var99;
        pardata.det1_pos2   = textdata.Var100;
        pardata.det1_pos3   = textdata.Var101;
       
        pardata.det2_pos1   = textdata.Var102;
        pardata.det2_pos2   = textdata.Var103;
        pardata.det2_pos3   = textdata.Var104;
        
        pardata.det3_pos1   = textdata.Var105;
        pardata.det3_pos2   = textdata.Var106;
        pardata.det3_pos3   = textdata.Var107;
        
        pardata.det4_pos1   = textdata.Var108;
        pardata.det4_pos2   = textdata.Var109;
        pardata.det4_pos3   = textdata.Var110;
        
        pardata.det5_pos1   = textdata.Var111;
        pardata.det5_pos2   = textdata.Var112;
        pardata.det5_pos3   = textdata.Var113;
        
        pardata.det6_pos1   = textdata.Var114;
        pardata.det6_pos2   = textdata.Var115;
        pardata.det6_pos3   = textdata.Var116;
      
        pardata.det7_pos1   = textdata.Var117;
        pardata.det7_pos2   = textdata.Var118;
        pardata.det7_pos3   = textdata.Var119;
       
        pardata.det8_pos1   = textdata.Var120;
        pardata.det8_pos2   = textdata.Var121;
        pardata.det8_pos3   = textdata.Var122;
       
        pardata.det9_pos1   = textdata.Var123;
        pardata.det9_pos2   = textdata.Var124;
        pardata.det9_pos3   = textdata.Var125;
        
        pardata.det10_pos1  = textdata.Var126;
        pardata.det10_pos2  = textdata.Var127;
        pardata.det10_pos3  = textdata.Var128;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var129;
        pardata.hex_pos2    = textdata.Var130;
        pardata.hex_pos3    = textdata.Var131;
        pardata.hex_pos4    = textdata.Var132;
        pardata.hex_pos5    = textdata.Var133;
        pardata.hex_pos6    = textdata.Var134;
        pardata.hex_pos7    = textdata.Var135;
       
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var136;
        pardata.slit1_V_pos     = textdata.Var137;
        pardata.slit1_H_size    = textdata.Var138;
        pardata.slit1_H_pos     = textdata.Var139;
       
        pardata.slit2_V_size    = textdata.Var140;
        pardata.slit2_V_pos     = textdata.Var141;
        pardata.slit2_H_size    = textdata.Var142;
        pardata.slit2_H_pos     = textdata.Var143;
        
        pardata.slit3_V_size    = textdata.Var144;
        pardata.slit3_V_pos     = textdata.Var145;
        pardata.slit3_H_size    = textdata.Var146;
        pardata.slit3_H_pos     = textdata.Var147;
        
        pardata.slit4_V_size    = textdata.Var148;
        pardata.slit4_V_pos     = textdata.Var149;
        pardata.slit4_H_size    = textdata.Var150;
        pardata.slit4_H_pos     = textdata.Var151;
        
        pardata.slit5_V_size    = textdata.Var152;
        pardata.slit5_V_pos     = textdata.Var153;
        pardata.slit5_H_size    = textdata.Var154;
        pardata.slit5_H_pos     = textdata.Var155;
       
        pardata.slit6_V_size    = textdata.Var156;
        pardata.slit6_V_pos     = textdata.Var157;
        pardata.slit6_H_size    = textdata.Var158;
        pardata.slit6_H_pos     = textdata.Var159;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var160;
        pardata.lens1_pos2  = textdata.Var161;
       
        pardata.lens2_pos1  = textdata.Var162;
        pardata.lens2_pos2  = textdata.Var163;
       
        pardata.lens3_pos1  = textdata.Var164;
        pardata.lens3_pos2  = textdata.Var165;
        
        pardata.lens4_pos1  = textdata.Var166;
        pardata.lens4_pos2  = textdata.Var167;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var168;
        pardata.encoder2    = textdata.Var169;
        pardata.encoder3    = textdata.Var170;
        pardata.encoder4    = textdata.Var171;
        pardata.encoder5    = textdata.Var172;
        pardata.encoder6    = textdata.Var173;
        pardata.encoder7    = textdata.Var174;
        pardata.encoder8    = textdata.Var175;
        pardata.encoder9    = textdata.Var176;
        pardata.encoder10   = textdata.Var177;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var178;
        pardata.ev2     = textdata.Var179;
        pardata.ev3     = textdata.Var180;
        pardata.ev4     = textdata.Var181;
        pardata.ev5     = textdata.Var182;
        pardata.ev6     = textdata.Var183;
        pardata.ev7     = textdata.Var184;
        pardata.ev8     = textdata.Var185;
        pardata.ev9     = textdata.Var186;
        pardata.ev10    = textdata.Var187;
    case 'ej_apr22'
        %%% READ IN DATA
        opts        = detectImportOptions(fname,'FileType','text', 'delimiter', ',');
        textdata    = readtable(fname, opts);

        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_frames_per_file    = textdata.Var42;
        pardata.det5_fnum               = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
        
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
        pardata.det10_fname             = textdata.Var61;
        pardata.det10_fnum              = textdata.Var62;
        pardata.det10_frames_per_file   = textdata.Var63;
        pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var65;
        pardata.scaler1_units   = textdata.Var66;
        pardata.scaler2_val     = textdata.Var67;
        pardata.scaler2_units   = textdata.Var68;
        pardata.scaler3_val     = textdata.Var69;
        pardata.scaler3_units   = textdata.Var70;
        pardata.scaler4_val     = textdata.Var71;
        pardata.scaler4_units   = textdata.Var72;
        pardata.scaler5_val     = textdata.Var73;
        pardata.scaler5_units   = textdata.Var74;
        pardata.scaler6_val     = textdata.Var75;
        pardata.scaler6_units   = textdata.Var76;
        pardata.scaler7_val     = textdata.Var77;
        pardata.scaler7_units   = textdata.Var78;
        pardata.scaler8_val     = textdata.Var79;
        pardata.scaler8_units   = textdata.Var80;
        pardata.scaler9_val     = textdata.Var81;
        pardata.scaler9_units   = textdata.Var82;
        pardata.scaler10_val    = textdata.Var83;
        pardata.scaler10_units  = textdata.Var84;
        pardata.scaler11_val    = textdata.Var85;
        pardata.scaler11_units  = textdata.Var86;
        pardata.scaler12_val    = textdata.Var87;
        pardata.scaler12_units  = textdata.Var88;
        
        pardata.samXE       = textdata.Var89;
        pardata.samYE       = textdata.Var90;
        pardata.samZE       = textdata.Var91;
        pardata.samRXE      = textdata.Var92;
        pardata.samRYE      = textdata.Var93;
        pardata.samRZE      = textdata.Var94;
        pardata.rotXE       = textdata.Var95;
        pardata.rotZE       = textdata.Var96;
        pardata.tblYE1      = textdata.Var97;
        pardata.tblYE2      = textdata.Var98;
        pardata.tblYE3      = textdata.Var99;
        pardata.mamRZE      = textdata.Var100;
        pardata.aerotblTiltX  = textdata.Var101;
        pardata.aerotblTiltX0 = textdata.Var102;
        pardata.aerotblTiltZ  = textdata.Var103;
        pardata.aerotblTiltZ0 = textdata.Var104;
        
        pardata.mtsXE       = textdata.Var105;
        pardata.mtsYE       = textdata.Var106;
        pardata.mtsZE       = textdata.Var107;
        pardata.mtsRXE      = textdata.Var108;
        pardata.mtsRYE      = textdata.Var109;
        pardata.mtsRZE      = textdata.Var110;
        pardata.mtsrotXE    = textdata.Var111;
        pardata.mtsrotZE    = textdata.Var112;
        pardata.mtstblYE1   = textdata.Var113;
        pardata.mtstblYE2   = textdata.Var114;
        pardata.mtstblYE3   = textdata.Var115;
        pardata.mtstblTiltX     = textdata.Var116;
        pardata.mtstblTiltX0    = textdata.Var117;
        pardata.mtstblTiltZ     = textdata.Var118;
        pardata.mtstblTiltZ0    = textdata.Var119;
        
        pardata.samXC       = textdata.Var120;
        pardata.samYC       = textdata.Var121;
        pardata.samZC       = textdata.Var122;
        pardata.samRXC      = textdata.Var123;
        pardata.samRYC      = textdata.Var124;
        pardata.samRZC      = textdata.Var125;
        pardata.rotXC       = textdata.Var126;
        pardata.rotZC       = textdata.Var127;
        pardata.tblYC1      = textdata.Var128;
        pardata.tblYC2      = textdata.Var129;
        pardata.tblYC3      = textdata.Var130;
        pardata.mamRZC      = textdata.Var131;
        
        pardata.samXrams3   = textdata.Var132;
        pardata.samYrams3   = textdata.Var133;
        pardata.samZrams3   = textdata.Var134;
        pardata.samRXrams3  = textdata.Var135;
        pardata.samRYrams3  = textdata.Var136;
        pardata.samRZrams3  = textdata.Var137;
        pardata.rotYTrams3  = textdata.Var138;
        pardata.rotYBrams3  = textdata.Var139;
        pardata.tblYrams3   = textdata.Var140;
        pardata.offsetrams3 = textdata.Var141;
        
        pardata.sampos_ev1  = textdata.Var142;
        pardata.sampos_ev2  = textdata.Var143;
        pardata.sampos_ev3  = textdata.Var144;
        pardata.sampos_ev4  = textdata.Var145;
        pardata.sampos_ev5  = textdata.Var146;
        pardata.sampos_ev6  = textdata.Var147;
        pardata.sampos_ev7  = textdata.Var148;
        pardata.sampos_ev8  = textdata.Var149;
        pardata.sampos_ev9  = textdata.Var150;
        pardata.sampos_ev10 = textdata.Var151;
        pardata.sampos_ev11 = textdata.Var152;
        pardata.sampos_ev12 = textdata.Var153;
        
        pardata.det1_pos1   = textdata.Var154;
        pardata.det1_pos2   = textdata.Var155;
        pardata.det1_pos3   = textdata.Var156;
       
        pardata.det2_pos1   = textdata.Var157;
        pardata.det2_pos2   = textdata.Var158;
        pardata.det2_pos3   = textdata.Var159;
        
        pardata.det3_pos1   = textdata.Var160;
        pardata.det3_pos2   = textdata.Var161;
        pardata.det3_pos3   = textdata.Var162;
        
        pardata.det4_pos1   = textdata.Var163;
        pardata.det4_pos2   = textdata.Var164;
        pardata.det4_pos3   = textdata.Var165;
        
        pardata.det5_pos1   = textdata.Var166;
        pardata.det5_pos2   = textdata.Var167;
        pardata.det5_pos3   = textdata.Var168;
        
        pardata.det6_pos1   = textdata.Var169;
        pardata.det6_pos2   = textdata.Var170;
        pardata.det6_pos3   = textdata.Var171;
      
        pardata.det7_pos1   = textdata.Var172;
        pardata.det7_pos2   = textdata.Var173;
        pardata.det7_pos3   = textdata.Var174;
       
        pardata.det8_pos1   = textdata.Var175;
        pardata.det8_pos2   = textdata.Var176;
        pardata.det8_pos3   = textdata.Var177;
       
        pardata.det9_pos1   = textdata.Var178;
        pardata.det9_pos2   = textdata.Var179;
        pardata.det9_pos3   = textdata.Var180;
        
        pardata.det10_pos1  = textdata.Var181;
        pardata.det10_pos2  = textdata.Var182;
        pardata.det10_pos3  = textdata.Var183;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var184;
        pardata.hex_pos2    = textdata.Var185;
        pardata.hex_pos3    = textdata.Var186;
        pardata.hex_pos4    = textdata.Var187;
        pardata.hex_pos5    = textdata.Var188;
        pardata.hex_pos6    = textdata.Var189;
        pardata.hex_pos7    = textdata.Var190;
       
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var191;
        pardata.slit1_V_pos     = textdata.Var192;
        pardata.slit1_H_size    = textdata.Var193;
        pardata.slit1_H_pos     = textdata.Var194;
       
        pardata.slit2_V_size    = textdata.Var195;
        pardata.slit2_V_pos     = textdata.Var196;
        pardata.slit2_H_size    = textdata.Var197;
        pardata.slit2_H_pos     = textdata.Var198;
        
        pardata.slit3_V_size    = textdata.Var199;
        pardata.slit3_V_pos     = textdata.Var200;
        pardata.slit3_H_size    = textdata.Var201;
        pardata.slit3_H_pos     = textdata.Var202;
        
        pardata.slit4_V_size    = textdata.Var203;
        pardata.slit4_V_pos     = textdata.Var204;
        pardata.slit4_H_size    = textdata.Var205;
        pardata.slit4_H_pos     = textdata.Var206;
        
        pardata.slit5_V_size    = textdata.Var207;
        pardata.slit5_V_pos     = textdata.Var208;
        pardata.slit5_H_size    = textdata.Var209;
        pardata.slit5_H_pos     = textdata.Var210;
       
        pardata.slit6_V_size    = textdata.Var211;
        pardata.slit6_V_pos     = textdata.Var212;
        pardata.slit6_H_size    = textdata.Var213;
        pardata.slit6_H_pos     = textdata.Var214;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var215;
        pardata.lens1_pos2  = textdata.Var216;
        
        pardata.lens2_pos1  = textdata.Var217;
        pardata.lens2_pos2  = textdata.Var218;
       
        pardata.lens3_pos1  = textdata.Var219;
        pardata.lens3_pos2  = textdata.Var220;
        
        pardata.lens4_pos1  = textdata.Var221;
        pardata.lens4_pos2  = textdata.Var222;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var223;
        pardata.encoder2    = textdata.Var224;
        pardata.encoder3    = textdata.Var225;
        pardata.encoder4    = textdata.Var226;
        pardata.encoder5    = textdata.Var227;
        pardata.encoder6    = textdata.Var228;
        pardata.encoder7    = textdata.Var229;
        pardata.encoder8    = textdata.Var230;
        pardata.encoder9    = textdata.Var231;
        pardata.encoder10   = textdata.Var232;
        pardata.encoder11   = textdata.Var233;
        pardata.encoder12   = textdata.Var234;
        pardata.encoder13   = textdata.Var235;
        pardata.encoder14   = textdata.Var236;
        pardata.encoder15   = textdata.Var237;
        pardata.encoder16   = textdata.Var238;
        pardata.encoder17   = textdata.Var239;
        pardata.encoder18   = textdata.Var240;
        pardata.encoder19   = textdata.Var241;
        pardata.encoder20   = textdata.Var242;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var243;
        pardata.ev2     = textdata.Var244;
        pardata.ev3     = textdata.Var245;
        pardata.ev4     = textdata.Var246;
        pardata.ev5     = textdata.Var247;
        pardata.ev6     = textdata.Var248;
        pardata.ev7     = textdata.Var249;
        pardata.ev8     = textdata.Var250;
        pardata.ev9     = textdata.Var251;
        pardata.ev10    = textdata.Var252;
    case 'saxs_waxs_fmt_fastpar_v6'
        %%% READ IN DATA
        opts        = detectImportOptions(fname, 'FileType', 'text');
        textdata    = readtable(fname, opts);
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_B_pos = textdata.Var13;
        pardata.atten_C_pos = textdata.Var14;
        
        pardata.det_type        = textdata.Var15;
        pardata.scan_mtr        = textdata.Var16;
        pardata.scan_ini        = textdata.Var17;
        pardata.scan_fin        = textdata.Var18;
        
        pardata.slave_mtr        = textdata.Var19;
        pardata.slave_mtr_ini    = textdata.Var20;
        pardata.slave_mtr_fin    = textdata.Var21;
        
        pardata.scan_nframes    = textdata.Var22;
        pardata.imgprefix       = textdata.Var23;
        pardata.imnum_ini       = textdata.Var24;
        pardata.imnum_fin       = textdata.Var25;
        
        pardata.det1_fname              = textdata.Var26;
        pardata.det1_fnum               = textdata.Var27;
        pardata.det1_frames_per_file    = textdata.Var28;
        pardata.det1_time_per_frame     = textdata.Var29;
        
        pardata.det2_fname              = textdata.Var30;
        pardata.det2_fnum               = textdata.Var31;
        pardata.det2_frames_per_file    = textdata.Var32;
        pardata.det2_time_per_frame     = textdata.Var33;
        
        pardata.det3_fname              = textdata.Var34;
        pardata.det3_fnum               = textdata.Var35;
        pardata.det3_frames_per_file    = textdata.Var36;
        pardata.det3_time_per_frame     = textdata.Var37;
       
        pardata.det4_fname              = textdata.Var38;
        pardata.det4_fnum               = textdata.Var39;
        pardata.det4_frames_per_file    = textdata.Var40;
        pardata.det4_time_per_frame     = textdata.Var41;
       
        pardata.det5_fname              = textdata.Var42;
        pardata.det5_frames_per_file    = textdata.Var43;
        pardata.det5_fnum               = textdata.Var44;
        pardata.det5_time_per_frame     = textdata.Var45;
      
        pardata.det6_fname              = textdata.Var46;
        pardata.det6_fnum               = textdata.Var47;
        pardata.det6_frames_per_file    = textdata.Var48;
        pardata.det6_time_per_frame     = textdata.Var49;
       
        pardata.det7_fname              = textdata.Var50;
        pardata.det7_fnum               = textdata.Var51;
        pardata.det7_frames_per_file    = textdata.Var52;
        pardata.det7_time_per_frame     = textdata.Var53;
       
        pardata.det8_fname              = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_fnum               = textdata.Var56;
        pardata.det8_time_per_frame     = textdata.Var57;
       
        pardata.det9_fname              = textdata.Var58;
        pardata.det9_fnum               = textdata.Var59;
        pardata.det9_frames_per_file    = textdata.Var60;
        pardata.det9_time_per_frame     = textdata.Var61;
        
        pardata.det10_fname             = textdata.Var62;
        pardata.det10_fnum              = textdata.Var63;
        pardata.det10_frames_per_file   = textdata.Var64;
        pardata.det10_time_per_frame    = textdata.Var65;
       
        pardata.scaler1_val     = textdata.Var66;
        pardata.scaler1_units   = textdata.Var67;
        pardata.scaler2_val     = textdata.Var68;
        pardata.scaler2_units   = textdata.Var69;
        pardata.scaler3_val     = textdata.Var70;
        pardata.scaler3_units   = textdata.Var71;
        pardata.scaler4_val     = textdata.Var72;
        pardata.scaler4_units   = textdata.Var73;
        pardata.scaler5_val     = textdata.Var74;
        pardata.scaler5_units   = textdata.Var75;
        pardata.scaler6_val     = textdata.Var76;
        pardata.scaler6_units   = textdata.Var77;
        pardata.scaler7_val     = textdata.Var78;
        pardata.scaler7_units   = textdata.Var79;
        pardata.scaler8_val     = textdata.Var80;
        pardata.scaler8_units   = textdata.Var81;
        pardata.scaler9_val     = textdata.Var82;
        pardata.scaler9_units   = textdata.Var83;
        pardata.scaler10_val    = textdata.Var84;
        pardata.scaler10_units  = textdata.Var85;
        pardata.scaler11_val    = textdata.Var86;
        pardata.scaler11_units  = textdata.Var87;
        pardata.scaler12_val    = textdata.Var88;
        pardata.scaler12_units  = textdata.Var89;
        
        pardata.samXE       = textdata.Var90;
        pardata.samYE       = textdata.Var91;
        pardata.samZE       = textdata.Var92;
        pardata.samRXE      = textdata.Var93;
        pardata.samRYE      = textdata.Var94;
        pardata.samRZE      = textdata.Var95;
        pardata.rotXE       = textdata.Var96;
        pardata.rotZE       = textdata.Var97;
        pardata.tblYE1      = textdata.Var98;
        pardata.tblYE2      = textdata.Var99;
        pardata.tblYE3      = textdata.Var100;
        pardata.mamRZE      = textdata.Var101;
        pardata.aerotblTiltX  = textdata.Var102;
        pardata.aerotblTiltX0 = textdata.Var103;
        pardata.aerotblTiltZ  = textdata.Var104;
        pardata.aerotblTiltZ0 = textdata.Var105;
       
        pardata.mtsXE       = textdata.Var106;
        pardata.mtsYE       = textdata.Var107;
        pardata.mtsZE       = textdata.Var108;
        pardata.mtsRXE      = textdata.Var109;
        pardata.mtsRYE      = textdata.Var110;
        pardata.mtsRZE      = textdata.Var111;
        pardata.mtsrotXE    = textdata.Var112;
        pardata.mtsrotZE    = textdata.Var113;
        pardata.mtstblYE1   = textdata.Var114;
        pardata.mtstblYE2   = textdata.Var115;
        pardata.mtstblYE3   = textdata.Var116;
        pardata.mtstblTiltX     = textdata.Var117;
        pardata.mtstblTiltX0    = textdata.Var118;
        pardata.mtstblTiltZ     = textdata.Var119;
        pardata.mtstblTiltZ0    = textdata.Var120;
        
        pardata.samXC       = textdata.Var121;
        pardata.samYC       = textdata.Var122;
        pardata.samZC       = textdata.Var123;
        pardata.samRXC      = textdata.Var124;
        pardata.samRYC      = textdata.Var125;
        pardata.samRZC      = textdata.Var126;
        pardata.rotXC       = textdata.Var127;
        pardata.rotZC       = textdata.Var128;
        pardata.tblYC1      = textdata.Var129;
        pardata.tblYC2      = textdata.Var130;
        pardata.tblYC3      = textdata.Var131;
        pardata.mamRZC      = textdata.Var132;
       
        pardata.samXrams3   = textdata.Var133;
        pardata.samYrams3   = textdata.Var134;
        pardata.samZrams3   = textdata.Var135;
        pardata.samRXrams3  = textdata.Var136;
        pardata.samRYrams3  = textdata.Var137;
        pardata.samRZrams3  = textdata.Var138;
        pardata.rotYTrams3  = textdata.Var139;
        pardata.rotYBrams3  = textdata.Var140;
        pardata.tblYrams3   = textdata.Var141;
        pardata.offsetrams3 = textdata.Var142;
      
        pardata.sampos_ev1  = textdata.Var143;
        pardata.sampos_ev2  = textdata.Var144;
        pardata.sampos_ev3  = textdata.Var145;
        pardata.sampos_ev4  = textdata.Var146;
        pardata.sampos_ev5  = textdata.Var147;
        pardata.sampos_ev6  = textdata.Var148;
        pardata.sampos_ev7  = textdata.Var149;
        pardata.sampos_ev8  = textdata.Var150;
        pardata.sampos_ev9  = textdata.Var151;
        pardata.sampos_ev10 = textdata.Var152;
        pardata.sampos_ev11 = textdata.Var153;
        pardata.sampos_ev12 = textdata.Var154;
       
        pardata.det1_pos1   = textdata.Var155;
        pardata.det1_pos2   = textdata.Var156;
        pardata.det1_pos3   = textdata.Var157;
       
        pardata.det2_pos1   = textdata.Var158;
        pardata.det2_pos2   = textdata.Var159;
        pardata.det2_pos3   = textdata.Var160;
       
        pardata.det3_pos1   = textdata.Var161;
        pardata.det3_pos2   = textdata.Var162;
        pardata.det3_pos3   = textdata.Var163;
        
        pardata.det4_pos1   = textdata.Var164;
        pardata.det4_pos2   = textdata.Var165;
        pardata.det4_pos3   = textdata.Var166;
        
        pardata.det5_pos1   = textdata.Var167;
        pardata.det5_pos2   = textdata.Var168;
        pardata.det5_pos3   = textdata.Var169;
        
        pardata.det6_pos1   = textdata.Var170;
        pardata.det6_pos2   = textdata.Var171;
        pardata.det6_pos3   = textdata.Var172;
      
        pardata.det7_pos1   = textdata.Var173;
        pardata.det7_pos2   = textdata.Var174;
        pardata.det7_pos3   = textdata.Var175;
       
        pardata.det8_pos1   = textdata.Var176;
        pardata.det8_pos2   = textdata.Var177;
        pardata.det8_pos3   = textdata.Var178;
       
        pardata.det9_pos1   = textdata.Var179;
        pardata.det9_pos2   = textdata.Var180;
        pardata.det9_pos3   = textdata.Var181;
        
        pardata.det10_pos1  = textdata.Var182;
        pardata.det10_pos2  = textdata.Var183;
        pardata.det10_pos3  = textdata.Var184;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var185;
        pardata.hex_pos2    = textdata.Var186;
        pardata.hex_pos3    = textdata.Var187;
        pardata.hex_pos4    = textdata.Var188;
        pardata.hex_pos5    = textdata.Var189;
        pardata.hex_pos6    = textdata.Var190;
        pardata.hex_pos7    = textdata.Var191;
      
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var192;
        pardata.slit1_V_pos     = textdata.Var193;
        pardata.slit1_H_size    = textdata.Var194;
        pardata.slit1_H_pos     = textdata.Var195;
      
        pardata.slit2_V_size    = textdata.Var196;
        pardata.slit2_V_pos     = textdata.Var197;
        pardata.slit2_H_size    = textdata.Var198;
        pardata.slit2_H_pos     = textdata.Var199;
        
        pardata.slit3_V_size    = textdata.Var200;
        pardata.slit3_V_pos     = textdata.Var201;
        pardata.slit3_H_size    = textdata.Var202;
        pardata.slit3_H_pos     = textdata.Var203;
       
        pardata.slit4_V_size    = textdata.Var204;
        pardata.slit4_V_pos     = textdata.Var205;
        pardata.slit4_H_size    = textdata.Var206;
        pardata.slit4_H_pos     = textdata.Var207;
       
        pardata.slit5_V_size    = textdata.Var208;
        pardata.slit5_V_pos     = textdata.Var209;
        pardata.slit5_H_size    = textdata.Var210;
        pardata.slit5_H_pos     = textdata.Var211;
       
        pardata.slit6_V_size    = textdata.Var212;
        pardata.slit6_V_pos     = textdata.Var213;
        pardata.slit6_H_size    = textdata.Var214;
        pardata.slit6_H_pos     = textdata.Var215;
       
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var216;
        pardata.lens1_pos2  = textdata.Var217;
       
        pardata.lens2_pos1  = textdata.Var218;
        pardata.lens2_pos2  = textdata.Var219;
       
        pardata.lens3_pos1  = textdata.Var220;
        pardata.lens3_pos2  = textdata.Var221;
        
        pardata.lens4_pos1  = textdata.Var222;
        pardata.lens4_pos2  = textdata.Var223;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var224;
        pardata.encoder2    = textdata.Var225;
        pardata.encoder3    = textdata.Var226;
        pardata.encoder4    = textdata.Var227;
        pardata.encoder5    = textdata.Var228;
        pardata.encoder6    = textdata.Var229;
        pardata.encoder7    = textdata.Var230;
        pardata.encoder8    = textdata.Var231;
        pardata.encoder9    = textdata.Var232;
        pardata.encoder10   = textdata.Var233;
        pardata.encoder11   = textdata.Var234;
        pardata.encoder12   = textdata.Var235;
        pardata.encoder13   = textdata.Var236;
        pardata.encoder14   = textdata.Var237;
        pardata.encoder15   = textdata.Var238;
        pardata.encoder16   = textdata.Var239;
        pardata.encoder17   = textdata.Var240;
        pardata.encoder18   = textdata.Var241;
        pardata.encoder19   = textdata.Var242;
        pardata.encoder20   = textdata.Var243;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var244;
        pardata.ev2     = textdata.Var245;
        pardata.ev3     = textdata.Var246;
        pardata.ev4     = textdata.Var247;
        pardata.ev5     = textdata.Var248;
        pardata.ev6     = textdata.Var249;
        pardata.ev7     = textdata.Var250;
        pardata.ev8     = textdata.Var251;
        pardata.ev9     = textdata.Var252;
        pardata.ev10    = textdata.Var253;
    case 'saxs_waxs_fmt_fastpar_v5'
        %%% READ IN DATA
        opts        = detectImportOptions(fname, 'FileType', 'text');
        textdata    = readtable(fname, opts);
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_frames_per_file    = textdata.Var42;
        pardata.det5_fnum               = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
        
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
        pardata.det10_fname             = textdata.Var61;
        pardata.det10_fnum              = textdata.Var62;
        pardata.det10_frames_per_file   = textdata.Var63;
        pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var65;
        pardata.scaler1_units   = textdata.Var66;
        pardata.scaler2_val     = textdata.Var67;
        pardata.scaler2_units   = textdata.Var68;
        pardata.scaler3_val     = textdata.Var69;
        pardata.scaler3_units   = textdata.Var70;
        pardata.scaler4_val     = textdata.Var71;
        pardata.scaler4_units   = textdata.Var72;
        pardata.scaler5_val     = textdata.Var73;
        pardata.scaler5_units   = textdata.Var74;
        pardata.scaler6_val     = textdata.Var75;
        pardata.scaler6_units   = textdata.Var76;
        pardata.scaler7_val     = textdata.Var77;
        pardata.scaler7_units   = textdata.Var78;
        pardata.scaler8_val     = textdata.Var79;
        pardata.scaler8_units   = textdata.Var80;
        pardata.scaler9_val     = textdata.Var81;
        pardata.scaler9_units   = textdata.Var82;
        pardata.scaler10_val    = textdata.Var83;
        pardata.scaler10_units  = textdata.Var84;
        pardata.scaler11_val    = textdata.Var85;
        pardata.scaler11_units  = textdata.Var86;
        pardata.scaler12_val    = textdata.Var87;
        pardata.scaler12_units  = textdata.Var88;
        
        pardata.samXE       = textdata.Var89;
        pardata.samYE       = textdata.Var90;
        pardata.samZE       = textdata.Var91;
        pardata.samRXE      = textdata.Var92;
        pardata.samRYE      = textdata.Var93;
        pardata.samRZE      = textdata.Var94;
        pardata.rotXE       = textdata.Var95;
        pardata.rotZE       = textdata.Var96;
        pardata.tblYE1      = textdata.Var97;
        pardata.tblYE2      = textdata.Var98;
        pardata.tblYE3      = textdata.Var99;
        pardata.mamRZE      = textdata.Var100;
        pardata.aerotblTiltX  = textdata.Var101;
        pardata.aerotblTiltX0 = textdata.Var102;
        pardata.aerotblTiltZ  = textdata.Var103;
        pardata.aerotblTiltZ0 = textdata.Var104;
        
        pardata.mtsXE       = textdata.Var105;
        pardata.mtsYE       = textdata.Var106;
        pardata.mtsZE       = textdata.Var107;
        pardata.mtsRXE      = textdata.Var108;
        pardata.mtsRYE      = textdata.Var109;
        pardata.mtsRZE      = textdata.Var110;
        pardata.mtsrotXE    = textdata.Var111;
        pardata.mtsrotZE    = textdata.Var112;
        pardata.mtstblYE1   = textdata.Var113;
        pardata.mtstblYE2   = textdata.Var114;
        pardata.mtstblYE3   = textdata.Var115;
        pardata.mtstblTiltX     = textdata.Var116;
        pardata.mtstblTiltX0    = textdata.Var117;
        pardata.mtstblTiltZ     = textdata.Var118;
        pardata.mtstblTiltZ0    = textdata.Var119;
        
        pardata.samXC       = textdata.Var120;
        pardata.samYC       = textdata.Var121;
        pardata.samZC       = textdata.Var122;
        pardata.samRXC      = textdata.Var123;
        pardata.samRYC      = textdata.Var124;
        pardata.samRZC      = textdata.Var125;
        pardata.rotXC       = textdata.Var126;
        pardata.rotZC       = textdata.Var127;
        pardata.tblYC1      = textdata.Var128;
        pardata.tblYC2      = textdata.Var129;
        pardata.tblYC3      = textdata.Var130;
        pardata.mamRZC      = textdata.Var131;
        
        pardata.samXrams3   = textdata.Var132;
        pardata.samYrams3   = textdata.Var133;
        pardata.samZrams3   = textdata.Var134;
        pardata.samRXrams3  = textdata.Var135;
        pardata.samRYrams3  = textdata.Var136;
        pardata.samRZrams3  = textdata.Var137;
        pardata.rotYTrams3  = textdata.Var138;
        pardata.rotYBrams3  = textdata.Var139;
        pardata.tblYrams3   = textdata.Var140;
        pardata.offsetrams3 = textdata.Var141;
        
        pardata.sampos_ev1  = textdata.Var142;
        pardata.sampos_ev2  = textdata.Var143;
        pardata.sampos_ev3  = textdata.Var144;
        pardata.sampos_ev4  = textdata.Var145;
        pardata.sampos_ev5  = textdata.Var146;
        pardata.sampos_ev6  = textdata.Var147;
        pardata.sampos_ev7  = textdata.Var148;
        pardata.sampos_ev8  = textdata.Var149;
        pardata.sampos_ev9  = textdata.Var150;
        pardata.sampos_ev10 = textdata.Var151;
        pardata.sampos_ev11 = textdata.Var152;
        pardata.sampos_ev12 = textdata.Var153;
        
        pardata.det1_pos1   = textdata.Var154;
        pardata.det1_pos2   = textdata.Var155;
        pardata.det1_pos3   = textdata.Var156;
       
        pardata.det2_pos1   = textdata.Var157;
        pardata.det2_pos2   = textdata.Var158;
        pardata.det2_pos3   = textdata.Var159;
        
        pardata.det3_pos1   = textdata.Var160;
        pardata.det3_pos2   = textdata.Var161;
        pardata.det3_pos3   = textdata.Var162;
        
        pardata.det4_pos1   = textdata.Var163;
        pardata.det4_pos2   = textdata.Var164;
        pardata.det4_pos3   = textdata.Var165;
        
        pardata.det5_pos1   = textdata.Var166;
        pardata.det5_pos2   = textdata.Var167;
        pardata.det5_pos3   = textdata.Var168;
        
        pardata.det6_pos1   = textdata.Var169;
        pardata.det6_pos2   = textdata.Var170;
        pardata.det6_pos3   = textdata.Var171;
      
        pardata.det7_pos1   = textdata.Var172;
        pardata.det7_pos2   = textdata.Var173;
        pardata.det7_pos3   = textdata.Var174;
       
        pardata.det8_pos1   = textdata.Var175;
        pardata.det8_pos2   = textdata.Var176;
        pardata.det8_pos3   = textdata.Var177;
       
        pardata.det9_pos1   = textdata.Var178;
        pardata.det9_pos2   = textdata.Var179;
        pardata.det9_pos3   = textdata.Var180;
        
        pardata.det10_pos1  = textdata.Var181;
        pardata.det10_pos2  = textdata.Var182;
        pardata.det10_pos3  = textdata.Var183;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var184;
        pardata.hex_pos2    = textdata.Var185;
        pardata.hex_pos3    = textdata.Var186;
        pardata.hex_pos4    = textdata.Var187;
        pardata.hex_pos5    = textdata.Var188;
        pardata.hex_pos6    = textdata.Var189;
        pardata.hex_pos7    = textdata.Var190;
       
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var191;
        pardata.slit1_V_pos     = textdata.Var192;
        pardata.slit1_H_size    = textdata.Var193;
        pardata.slit1_H_pos     = textdata.Var194;
       
        pardata.slit2_V_size    = textdata.Var195;
        pardata.slit2_V_pos     = textdata.Var196;
        pardata.slit2_H_size    = textdata.Var197;
        pardata.slit2_H_pos     = textdata.Var198;
        
        pardata.slit3_V_size    = textdata.Var199;
        pardata.slit3_V_pos     = textdata.Var200;
        pardata.slit3_H_size    = textdata.Var201;
        pardata.slit3_H_pos     = textdata.Var202;
        
        pardata.slit4_V_size    = textdata.Var203;
        pardata.slit4_V_pos     = textdata.Var204;
        pardata.slit4_H_size    = textdata.Var205;
        pardata.slit4_H_pos     = textdata.Var206;
        
        pardata.slit5_V_size    = textdata.Var207;
        pardata.slit5_V_pos     = textdata.Var208;
        pardata.slit5_H_size    = textdata.Var209;
        pardata.slit5_H_pos     = textdata.Var210;
       
        pardata.slit6_V_size    = textdata.Var211;
        pardata.slit6_V_pos     = textdata.Var212;
        pardata.slit6_H_size    = textdata.Var213;
        pardata.slit6_H_pos     = textdata.Var214;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var215;
        pardata.lens1_pos2  = textdata.Var216;
        
        pardata.lens2_pos1  = textdata.Var217;
        pardata.lens2_pos2  = textdata.Var218;
       
        pardata.lens3_pos1  = textdata.Var219;
        pardata.lens3_pos2  = textdata.Var220;
        
        pardata.lens4_pos1  = textdata.Var221;
        pardata.lens4_pos2  = textdata.Var222;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var223;
        pardata.encoder2    = textdata.Var224;
        pardata.encoder3    = textdata.Var225;
        pardata.encoder4    = textdata.Var226;
        pardata.encoder5    = textdata.Var227;
        pardata.encoder6    = textdata.Var228;
        pardata.encoder7    = textdata.Var229;
        pardata.encoder8    = textdata.Var230;
        pardata.encoder9    = textdata.Var231;
        pardata.encoder10   = textdata.Var232;
        pardata.encoder11   = textdata.Var233;
        pardata.encoder12   = textdata.Var234;
        pardata.encoder13   = textdata.Var235;
        pardata.encoder14   = textdata.Var236;
        pardata.encoder15   = textdata.Var237;
        pardata.encoder16   = textdata.Var238;
        pardata.encoder17   = textdata.Var239;
        pardata.encoder18   = textdata.Var240;
        pardata.encoder19   = textdata.Var241;
        pardata.encoder20   = textdata.Var242;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var243;
        pardata.ev2     = textdata.Var244;
        pardata.ev3     = textdata.Var245;
        pardata.ev4     = textdata.Var246;
        pardata.ev5     = textdata.Var247;
        pardata.ev6     = textdata.Var248;
        pardata.ev7     = textdata.Var249;
        pardata.ev8     = textdata.Var250;
        pardata.ev9     = textdata.Var251;
        pardata.ev10    = textdata.Var252;
    case 'saxs_waxs_fmt_fastpar_wzhang_sep21'
        %%% STARTING brown_jun21 (after edits)
        %%% READ IN DATA
        opts        = detectImportOptions(fname,'FileType','text');
        textdata    = readtable(fname, opts);
        % keyboard
        % return
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.epoch_time  = textdata.Var6;
        pardata.integ_time  = textdata.Var7;
        pardata.Iring       = textdata.Var8;
        pardata.und_gap     = textdata.Var9;
        pardata.energy      = textdata.Var10;
        pardata.energy_cal  = textdata.Var11;
        pardata.foil_pos    = textdata.Var12;
        pardata.atten_pos   = textdata.Var13;
        
        pardata.det_type        = textdata.Var14;
        pardata.scan_mtr        = textdata.Var15;
        pardata.scan_ini        = textdata.Var16;
        pardata.scan_fin        = textdata.Var17;
        
        pardata.slave_mtr        = textdata.Var18;
        pardata.slave_mtr_ini    = textdata.Var19;
        pardata.slave_mtr_fin    = textdata.Var20;
        
        pardata.scan_nframes    = textdata.Var21;
        pardata.imgprefix       = textdata.Var22;
        pardata.imnum_ini       = textdata.Var23;
        pardata.imnum_fin       = textdata.Var24;
        
        pardata.det1_fname              = textdata.Var25;
        pardata.det1_fnum               = textdata.Var26;
        pardata.det1_frames_per_file    = textdata.Var27;
        pardata.det1_time_per_frame     = textdata.Var28;
        
        pardata.det2_fname              = textdata.Var29;
        pardata.det2_fnum               = textdata.Var30;
        pardata.det2_frames_per_file    = textdata.Var31;
        pardata.det2_time_per_frame     = textdata.Var32;
        
        pardata.det3_fname              = textdata.Var33;
        pardata.det3_fnum               = textdata.Var34;
        pardata.det3_frames_per_file    = textdata.Var35;
        pardata.det3_time_per_frame     = textdata.Var36;
       
        pardata.det4_fname              = textdata.Var37;
        pardata.det4_fnum               = textdata.Var38;
        pardata.det4_frames_per_file    = textdata.Var39;
        pardata.det4_time_per_frame     = textdata.Var40;
        
        pardata.det5_fname              = textdata.Var41;
        pardata.det5_frames_per_file    = textdata.Var42;
        pardata.det5_fnum               = textdata.Var43;
        pardata.det5_time_per_frame     = textdata.Var44;
       
        pardata.det6_fname              = textdata.Var45;
        pardata.det6_fnum               = textdata.Var46;
        pardata.det6_frames_per_file    = textdata.Var47;
        pardata.det6_time_per_frame     = textdata.Var48;
       
        pardata.det7_fname              = textdata.Var49;
        pardata.det7_fnum               = textdata.Var50;
        pardata.det7_frames_per_file    = textdata.Var51;
        pardata.det7_time_per_frame     = textdata.Var52;
        
        pardata.det8_fname              = textdata.Var53;
        pardata.det8_fnum               = textdata.Var54;
        pardata.det8_frames_per_file    = textdata.Var55;
        pardata.det8_time_per_frame     = textdata.Var56;
        
        pardata.det9_fname              = textdata.Var57;
        pardata.det9_fnum               = textdata.Var58;
        pardata.det9_frames_per_file    = textdata.Var59;
        pardata.det9_time_per_frame     = textdata.Var60;
        
%         pardata.det10_fname             = textdata.Var61;
%         pardata.det10_fnum              = textdata.Var62;
%         pardata.det10_frames_per_file   = textdata.Var63;
%         pardata.det10_time_per_frame    = textdata.Var64;
        
        pardata.scaler1_val     = textdata.Var63;
        pardata.scaler1_units   = textdata.Var64;
        pardata.scaler2_val     = textdata.Var65;
        pardata.scaler2_units   = textdata.Var66;
        pardata.scaler3_val     = textdata.Var67;
        pardata.scaler3_units   = textdata.Var68;
        pardata.scaler4_val     = textdata.Var69;
        pardata.scaler4_units   = textdata.Var70;
        pardata.scaler5_val     = textdata.Var71;
        pardata.scaler5_units   = textdata.Var72;
        pardata.scaler6_val     = textdata.Var73;
        pardata.scaler6_units   = textdata.Var74;
        pardata.scaler7_val     = textdata.Var75;
        pardata.scaler7_units   = textdata.Var76;
        pardata.scaler8_val     = textdata.Var77;
        pardata.scaler8_units   = textdata.Var78;
        pardata.scaler9_val     = textdata.Var79;
        pardata.scaler9_units   = textdata.Var80;
        pardata.scaler10_val    = textdata.Var81;
        pardata.scaler10_units  = textdata.Var82;
        pardata.scaler11_val    = textdata.Var83;
        pardata.scaler11_units  = textdata.Var84;
        pardata.scaler12_val    = textdata.Var85;
        pardata.scaler12_units  = textdata.Var86;
        
        pardata.samX        = textdata.Var87;
        pardata.samY        = textdata.Var88;
        pardata.samZ        = textdata.Var89;
        pardata.aX          = textdata.Var90;
        pardata.aY          = textdata.Var91;
        pardata.aZ          = textdata.Var92;
        pardata.samX2       = textdata.Var93;
        pardata.samY2       = textdata.Var94;
        pardata.samZ2       = textdata.Var95;
        pardata.samOther    = textdata.Var96;
        
        pardata.det1_pos1   = textdata.Var97;
        pardata.det1_pos2   = textdata.Var98;
        pardata.det1_pos3   = textdata.Var99;
       
        pardata.det2_pos1   = textdata.Var100;
        pardata.det2_pos2   = textdata.Var101;
        pardata.det2_pos3   = textdata.Var102;
        
        pardata.det3_pos1   = textdata.Var103;
        pardata.det3_pos2   = textdata.Var104;
        pardata.det3_pos3   = textdata.Var105;
        
        pardata.det4_pos1   = textdata.Var106;
        pardata.det4_pos2   = textdata.Var107;
        pardata.det4_pos3   = textdata.Var108;
        
        pardata.det5_pos1   = textdata.Var109;
        pardata.det5_pos2   = textdata.Var110;
        pardata.det5_pos3   = textdata.Var111;
        
        pardata.det6_pos1   = textdata.Var112;
        pardata.det6_pos2   = textdata.Var113;
        pardata.det6_pos3   = textdata.Var114;
      
        pardata.det7_pos1   = textdata.Var115;
        pardata.det7_pos2   = textdata.Var116;
        pardata.det7_pos3   = textdata.Var117;
       
        pardata.det8_pos1   = textdata.Var118;
        pardata.det8_pos2   = textdata.Var119;
        pardata.det8_pos3   = textdata.Var120;
       
        pardata.det9_pos1   = textdata.Var121;
        pardata.det9_pos2   = textdata.Var122;
        pardata.det9_pos3   = textdata.Var123;
        
        pardata.det10_pos1  = textdata.Var124;
        pardata.det10_pos2  = textdata.Var125;
        pardata.det10_pos3  = textdata.Var126;
        
        %%% LINE 8
        pardata.hex_pos1    = textdata.Var127;
        pardata.hex_pos2    = textdata.Var128;
        pardata.hex_pos3    = textdata.Var129;
        pardata.hex_pos4    = textdata.Var130;
        pardata.hex_pos5    = textdata.Var131;
        pardata.hex_pos6    = textdata.Var132;
        pardata.hex_pos7    = textdata.Var133;
       
        %%% LINE 9
        pardata.slit1_V_size    = textdata.Var134;
        pardata.slit1_V_pos     = textdata.Var135;
        pardata.slit1_H_size    = textdata.Var136;
        pardata.slit1_H_pos     = textdata.Var137;
       
        pardata.slit2_V_size    = textdata.Var138;
        pardata.slit2_V_pos     = textdata.Var139;
        pardata.slit2_H_size    = textdata.Var140;
        pardata.slit2_H_pos     = textdata.Var141;
        
        pardata.slit3_V_size    = textdata.Var142;
        pardata.slit3_V_pos     = textdata.Var143;
        pardata.slit3_H_size    = textdata.Var144;
        pardata.slit3_H_pos     = textdata.Var145;
        
        pardata.slit4_V_size    = textdata.Var146;
        pardata.slit4_V_pos     = textdata.Var147;
        pardata.slit4_H_size    = textdata.Var148;
        pardata.slit4_H_pos     = textdata.Var149;
        
        pardata.slit5_V_size    = textdata.Var150;
        pardata.slit5_V_pos     = textdata.Var151;
        pardata.slit5_H_size    = textdata.Var152;
        pardata.slit5_H_pos     = textdata.Var153;
       
        pardata.slit6_V_size    = textdata.Var154;
        pardata.slit6_V_pos     = textdata.Var155;
        pardata.slit6_H_size    = textdata.Var156;
        pardata.slit6_H_pos     = textdata.Var157;
        
        %%% LINE 10
        pardata.lens1_pos1  = textdata.Var158;
        pardata.lens1_pos2  = textdata.Var159;
       
        pardata.lens2_pos1  = textdata.Var160;
        pardata.lens2_pos2  = textdata.Var161;
       
        pardata.lens3_pos1  = textdata.Var162;
        pardata.lens3_pos2  = textdata.Var163;
        
        pardata.lens4_pos1  = textdata.Var164;
        pardata.lens4_pos2  = textdata.Var165;
       
        %%% LINE 11
        pardata.encoder1    = textdata.Var166;
        pardata.encoder2    = textdata.Var167;
        pardata.encoder3    = textdata.Var168;
        pardata.encoder4    = textdata.Var169;
        pardata.encoder5    = textdata.Var170;
        pardata.encoder6    = textdata.Var171;
        pardata.encoder7    = textdata.Var172;
        pardata.encoder8    = textdata.Var173;
        pardata.encoder9    = textdata.Var174;
        pardata.encoder10   = textdata.Var175;
        
        %%% LINE 12
        pardata.ev1     = textdata.Var176;
        pardata.ev2     = textdata.Var177;
        pardata.ev3     = textdata.Var178;
        pardata.ev4     = textdata.Var179;
        pardata.ev5     = textdata.Var180;
        pardata.ev6     = textdata.Var181;
        pardata.ev7     = textdata.Var182;
        pardata.ev8     = textdata.Var183;
        pardata.ev9     = textdata.Var184;
        pardata.ev10    = textdata.Var185;
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
    case 'mpe_ff_per_frame_v1'
        fmtstring   = [ ...
            '%s %s %d %s %d ' ...
            '%s %f %f %f %f %f %f %f %f %f ' ...
            '%f %f %f %f %f %f %f %f %f %f ' ...
            '%s %f %f %f %f ' ...
            '%s %d %d %d %f %f ' ...
            '%f %f %f %f'];

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
    case 'mpe_ff_per_frame_v2'
        %%% STARTING mpe1_oct20
        %%% READ IN DATA
        opts        = detectImportOptions(fname,'FileType','text');
        textdata    = readtable(fname, opts);
        
        %%% PARSE DATA
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.det_type    = textdata.Var6;
        pardata.det_fname   = textdata.Var7;
        
        pardata.scanmtr = textdata.Var9;
        pardata.scanini = textdata.Var10;
        pardata.scanfin = textdata.Var11;
        
        pardata.samX    = textdata.Var12;
        pardata.samY    = textdata.Var13;
        pardata.samZ    = textdata.Var14;
        
        pardata.aX      = textdata.Var15;
        pardata.aY      = textdata.Var17;
        pardata.aZ      = textdata.Var16;
        
        pardata.det1_fnum           = textdata.Var8;
        pardata.det1_time_per_frame = textdata.Var18;
        pardata.det1_frame_number   = textdata.Var19;
        
        pardata.scaler_count_time   = textdata.Var20;
        
        pardata.scaler1_val = textdata.Var21;
        pardata.scaler2_val = textdata.Var22;
        
        pardata.encoder1    = textdata.Var24;
        pardata.encoder2    = textdata.Var25;
        pardata.encoder3    = textdata.Var26;
        pardata.encoder4    = textdata.Var30;
        
        pardata.scaler3_val = textdata.Var27;
        pardata.scaler4_val = textdata.Var28;
        
        pardata.ev1 = textdata.Var29;
        
        pardata.samX2       = textdata.Var31;
        
        pardata.slit1_V_size    = textdata.Var33;
        pardata.slit1_V_pos     = textdata.Var34;
        pardata.slit1_H_size    = textdata.Var35;
        pardata.slit1_H_pos     = textdata.Var36;
        
        pardata.slit2_V_size    = textdata.Var37;
        pardata.slit2_V_pos     = textdata.Var38;
        pardata.slit2_H_size    = textdata.Var39;
        pardata.slit2_H_pos     = textdata.Var40;
        
        
        pardata.scanpos = textdata.Var12;
        
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
    case 'park_aug20'
        textdata    = readtable(fname, 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
         
        pardata.day     = textdata.Var1;
        pardata.month   = textdata.Var2;
        pardata.date    = textdata.Var3;
        pardata.time    = textdata.Var4;
        pardata.year    = textdata.Var5;
        
        pardata.det_type	= textdata.Var6;
        pardata.Iring       = textdata.Var7;
        
        pardata.aZ          = textdata.Var8;
        pardata.aY          = nan;
        pardata.aX          = textdata.Var9;
        
        pardata.samY        = textdata.Var10;
        pardata.samX        = textdata.Var11;
        pardata.samZ        = textdata.Var12;
        
        pardata.samX2       = textdata.Var13;
        pardata.samY2       = nan;
        pardata.samZ2       = nan;
        
        pardata.samOther    = nan;
        
        pardata.scanmtr     = textdata.Var16;
        pardata.scanini     = textdata.Var17;
        pardata.scanfin     = textdata.Var18;
        
        pardata.det1_frames_per_file    = textdata.Var19;
        pardata.det1_time_per_frame     = textdata.Var20;
        pardata.det1_fname              = textdata.Var21;
        pardata.det1_fnum               = textdata.Var22;
        pardata.integ_time              = textdata.Var20;
        
        pardata.scaler1_val     = textdata.Var24;
        pardata.scaler1_units   = nan;
        pardata.scaler2_val     = textdata.Var25;
        pardata.scaler2_units   = nan;
        pardata.scaler3_val     = textdata.Var26;
        pardata.scaler3_units   = nan;
        pardata.scaler4_val     = textdata.Var27;
        pardata.scaler4_units   = nan;
        pardata.scaler5_val     = textdata.Var28;
        pardata.scaler5_units   = nan;
        pardata.scaler6_val     = textdata.Var29;
        pardata.scaler6_units   = nan;
        pardata.scaler7_val     = textdata.Var30;
        pardata.scaler7_units   = nan;
        pardata.scaler8_val     = textdata.Var31;
        pardata.scaler8_units   = nan;
        pardata.scaler9_val     = nan;
        pardata.scaler9_units   = nan;
        pardata.scaler10_val    = nan;
        pardata.scaler10_units  = nan;
        
        pardata.encoder1    = textdata.Var32;
        pardata.encoder2    = textdata.Var33;
        pardata.encoder3    = textdata.Var34;
        pardata.encoder4    = textdata.Var35;
        pardata.encoder5    = textdata.Var36;
        pardata.encoder6    = textdata.Var37;
        pardata.encoder7    = textdata.Var38;
        pardata.encoder8    = textdata.Var39;
        pardata.encoder9    = textdata.Var40;
        pardata.encoder10   = textdata.Var41;
        
        pardata.ev1     = textdata.Var50;
        pardata.ev2     = textdata.Var51;
        pardata.ev3     = nan;
        pardata.ev4     = nan;
        pardata.ev5     = nan;
        pardata.ev6     = nan;
        pardata.ev7     = nan;
        pardata.ev8     = nan;
        pardata.ev9     = nan;
        pardata.ev10    = nan;
        
        pardata.slit1_V_size    = textdata.Var42;
        pardata.slit1_V_pos     = textdata.Var43;
        pardata.slit1_H_size    = textdata.Var44;
        pardata.slit1_H_pos     = textdata.Var45;
        
        pardata.slit2_V_size    = textdata.Var46;
        pardata.slit2_V_pos     = textdata.Var47;
        pardata.slit2_H_size    = textdata.Var48;
        pardata.slit2_H_pos     = textdata.Var49;
        
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
        
        pardata.foil_pos    = textdata.Var52;
        pardata.atten_pos   = textdata.Var53;
        pardata.atten_Cpos  = textdata.Var54;
        
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
    otherwise
        disp('format not implemented')
end
