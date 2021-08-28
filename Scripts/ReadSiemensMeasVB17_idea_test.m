clear all; close all;
% measfile = 'meas_MID43_MEGE_B0 _4ch_tuneup_FID17436.dat';
% measfile = 'meas_MID160_wip_fl_ute_FID7646.dat';
% measfile = '8channel_B1maps_phantom.dat';
% measfile = 'C:\MIDEA\N4_VB17A_LATEST_20090307\n4\pkg\MrServers\MrVista\Ice\MyIceFunctors\IceTemplate\MeasData\meas.dat';
% measfile = 'C:\IDEA\n4\pkg\MeasurementData\VB15\meas_MID11_gre_FID11691.dat';
% measfile = 'C:\IDEA\n4\pkg\MeasurementData\VB15\meas.dat';
% measfile = 'meas_MID160_wip_fl_ute_FID7646.dat';

% measfile = fullfile('Measdata', 'meas_MID43_MEGE_B0_4ch_tuneup_FID17436.dat');
% measfile = fullfile('Measdata', 'meas_MID160_wip_fl_ute_FID7646.dat');
% measfile = fullfile('Measdata', 'meas_MID154_wip_fl_ute_FID29428.dat');  
% measfile = fullfile('Measdata', '20110401', 'meas_MID75_wip_fl_ute_FID48512.dat');  

% measfile = fullfile('Measdata', 'meas_MID218_CV_Cine_prospectiv_11_segments_sense_4_FID19670.dat');
% measfile = fullfile('Measdata', 'meas_MID220_CV_Cine_retrospectiv_11_segments_sense_2_FID19672.dat');
% measfile = fullfile('Measdata', 'meas_MID221_TSE_turbo_factor_10_FID19673.dat');
% measfile = fullfile('Measdata', 'meas_MID225_TSE_shifted_FID19677.dat');
% measfile = fullfile('Measdata', 'meas_MID38_TSE_FID23798.dat');
% measfile = fullfile('Measdata', 'meas_MID39_tse_displace_FID23799.dat');
% measfile = fullfile('Measdata', 'meas_MID134_tse_displace_amp_faktor_3_FID24655.dat');
% measfile = fullfile('Measdata', 'meas_MID340_wip_fl_ute_19F.dat');
% measfile = fullfile('Measdata', 'displaced_TSE', 'meas_MID12_tse_displaced_singleshot_sense2_FID26405.dat');  
% measfile = fullfile('Measdata', 'displaced_TSE', 'meas_MID24_tse_cardio_displacedUFLARE_FID27543.dat');  
% measfile = fullfile('Measdata', 'displaced_TSE', 'meas_MID316_tse_cardio_split_echo_FID27931.dat');  
% measfile = fullfile('Measdata', 'displaced_TSE', 'meas_MID151_invivo_tse_displaced_multishot_sense2_FID26393.dat');  
% measfile = fullfile('Measdata', 'displaced_TSE', 'meas_MID150_invivo_tse_displaced_singleshot_sense2_FID26392.dat');  
% measfile = '/mridata1/mri_group/maolin_data/Yale/matlab/Measdata/magalie/meas_MID263_ep_seg_therm_ROgated_1sl_1dyn_FID58576.dat';  

measfile = '/mridata1/mri_group/maolin_data/Yale/matlab/Measdata/meas_MID250_FLASH_radial_me_64_FID35304.dat';
measfile = '/mridata1/mri_group/maolin_data/Yale/matlab/Measdata/meas.dat';

% [pa fa ex ve] = fileparts(pwd);
% [pa1 fa1 ex1 ve1] = fileparts(pa);
 
% measfile = fullfile(pa1, measfile);
 
% measfile = '/mridata1/mri_group/maolin_data/Temp/20110413/meas_MID46_wip_fl_ute_FID49448.dat';
% measfile = '/mridata1/mri_group/maolin_data/Temp/20110413/meas_MID50_T2_highres_FID49452.dat';
% measfile = '/mridata1/mri_group/maolin_data/Temp/20110413/meas_MID51_T2_highres_mSENSE_FID49453.dat';

% [rawdata loopcounters sMDH MrProt Yaps lc_names] = ReadSiemensMeasVB17(measfile, 'on');
[rawdata loopcounters sMDH lc_names] = ReadSiemensMeasVB17_idea(measfile, 'on');