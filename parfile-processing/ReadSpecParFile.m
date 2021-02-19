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
    
    
    
    case 'saxs_waxs_fmt_fastpar_v2'
        %%% STARTING park_aug20 (maybe even before when slave motors for fs was enabled)
        
        %%% READ IN DATA
        textdata    = readtable(fname, 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'split');
        
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
