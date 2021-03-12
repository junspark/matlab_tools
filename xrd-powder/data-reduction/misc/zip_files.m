clear all
close all
clc

froot{1}    = 'DP980_45_4_waxs_cont_loading';
froot{2}    = 'Ti_EOS_H_Ti64_900_sam1_waxs_cont_loading';
froot{3}    = 'DP980_90_5_waxs_cont_loading';
froot{4}    = 'Ti_EOS_H_Ti64_ab_sam1_waxs_cont_loading';
froot{5}    = 'LaB6_0pt3s';
froot{6}    = 'Ti_EOS_V_Ti44_400_sam1_waxs_cont_loading';
froot{7}    = 'LaB6_1s_200umx200um_beam';
froot{8}    = 'Ti_EOS_V_Ti44_ab_sam2_waxs_cont_loading';
froot{9}    = 'LaB6_1s_Mg';
froot{10}   = 'Ti_EOS_V_Ti64_400_sam2_waxs_cont_loading';
froot{11}    = 'MgAlZnCa_no2_sam1_waxs_cont_loading';
froot{12}    = 'Ti_EOS_V_Ti64_730_sam1_waxs_cont_loading';
froot{13}    = 'MgAlZnCa_no3_sam1_waxs_cont_loading';
froot{14}    = 'Ti_EOS_V_Ti64_900_sam1_waxs_cont_loading';
froot{15}    = 'MgAlZnCa_no5_sam1_waxs_cont_loading';
froot{16}    = 'Ti_EOS_V_Ti64_ab_sam2_waxs_cont_loading';
froot{17}    = 'MgAlZnCa_no9_sam1_waxs_cont_loading';
froot{18}    = 'Ti_HS_H_Ti64_400_sam1_waxs_cont_loading';
froot{19}    = 'MgAX138_sam1_waxs_cont_loading';
froot{20}    = 'Ti_HS_H_Ti64_ab_sam1_waxs_cont_loading';
froot{21}    = 'MgAX138_sam2_waxs_cont_loading';
froot{22}    = 'Ti_HS_V_Ti64_400_sam1_waxs_cont_loading';
froot{23}    = 'MgL3AX83_sam1_waxs_cont_loading';
froot{24}    = 'Ti_HS_V_Ti64_ab_sam1_waxs_cont_loading';
froot{25}    = 'static_Ti_EOS_V_Ti44_ab';
froot{26}    = 'Ti_MIT_H_waxs_cont_loading';
froot{27}    = 'static_Ti_EOS_V_Ti64_400';
froot{28}    = 'Ti_MIT_raw_waxs_cont_loading';
froot{29}    = 'static_Ti_EOS_V_Ti64_ab';
froot{30}    = 'Ti_EOS_H_Ti44_400_sam1_waxs_cont_loading';
froot{31}    = 'Ti_EOS_H_Ti44_ab_sam1_waxs_cont_loading';
froot{32}    = 'Ti_EOS_H_Ti64_400_sam2_waxs_cont_loading';
froot{33}    = 'Ti_EOS_H_Ti64_730_sam1_waxs_cont_loading';

pname_base  = '/home/beams/S1IDUSER/mnt/orthros/lywang_nov20_bc';
for iii = 1:1:length(froot)
    pname   = fullfile(pname_base, froot{iii});
    cmdstr  = sprintf('! zip %s %s/ge3/*.mat', pname, pname);
    eval(cmdstr)
end