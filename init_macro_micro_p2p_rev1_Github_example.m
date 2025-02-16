clear;
clc;
close all;
close all force;
app=NaN(1);  %%%%%%%%%This is to allow for Matlab Application integration.
format shortG
top_start_clock=clock;
folder1='C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\7GHz Non-Aggregate Ground';
cd(folder1)
addpath(folder1)
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Basic_Functions')
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\General_Movelist')  %%%%%%%%This is another Github repo
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\General_Terrestrial_Pathloss')  %%%%%%%%This is another Github repo
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Generic_Bugsplat')
addpath('C:\Local Matlab Data') %%%%%%%One Drive Error with mat files(These files are in the other github repos.
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'For the Macro/Micro and Fixed p2p non-aggregate coordination zones'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Base Station Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
load('aas_zero_elevation_data.mat','aas_zero_elevation_data')
toc;
%%%%1) Azimuth -180~~180
%%%2) Rural
%%%3) Suburban
%%%4) Urban
%%%%AAS Reduction in Gain to Max Gain (0dB is 0dB reduction, which equates to the make antenna gain of 25dB)
%%%%Need to normalize to zero after the "downtilt reductions" are calculated
%%%%To simplify the data, this is gain at the horizon. 50th Percentile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aas_zero_elevation_data(1,:)
%%%%%%%%%%Set all gains to 0dB, since there will be no azimuths
aas_zero_elevation_data(:,[2:4])=0;
aas_zero_elevation_data(1,:)
bs_down_tilt_reduction=abs(max(aas_zero_elevation_data(:,[2:4]))) %%%%%%%%Downtilt dB Value for Rural/Suburban/Urban
bs_down_tilt_reduction=0
norm_aas_zero_elevation_data=horzcat(aas_zero_elevation_data(:,1),aas_zero_elevation_data(:,[2:4])+bs_down_tilt_reduction);
max(norm_aas_zero_elevation_data(:,[2:4])) %%%%%This should be [0 0 0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%cell_p2p_data
%     % %1)Name,
%     % % 2)Lat
%     % %%3)Lon,
%     % % 4) Receiver Bandwidth (MHz),
%     % % 5) Antenna Height,
%     % % 6) Antenna Horizontal Beamwidth,
%     % % 7) Min Azimuth (Can be the same at Max Azimuth)
%     % % 8) Max Azimuth
%     % % 9) Antenna gain
%     % % 10) Rx NF
%     % % 11) I/N Ratio
%     % % 12) Main to side gain: 40dB
%     % % 13) FDR dB
%     % % 14) DPA Threshold [xx/100MHz]
%     % % 15) Required Pathloss

cell_p2p_data=cell(1,15);
cell_p2p_data{1}='Dodgers_Stadium';
cell_p2p_data{2}=34.082898;
cell_p2p_data{3}=-118.243378;
cell_p2p_data{4}=30; %%%%%%Mhz % P2P Receiver Bandwidth: [30 megaherz]
cell_p2p_data{5}=35;%%%%Meters % P2P Antenna Height: [35 meters]
cell_p2p_data{6}=2; % P2P Antenna 3dB Beamwidth: [2°]
cell_p2p_data{7}=135; %%%%%Min Azimuth
cell_p2p_data{8}=135; %%%%%Max Azimuth (This is the same because it's fixed)
cell_p2p_data{9}=40; % P2P Antenna Gain: [40dBi]
cell_p2p_data{10}=5; %%%%P2P Rx NF dB
cell_p2p_data{11}=-6; %%%%%I/N Ratio -6dB
cell_p2p_data{12}=45;%%%%%%%%Main to side gain: 45dB


cell_p2p_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Macrocell Example
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev=201; %%%%%%P2P Test: Github Example: MacroCell Dodgers
% grid_spacing=1;  %%%%km:
% bs_eirp=82%%%%horzcat(82,82,82); %%%%%EIRP [dBm/100MHz] for Rural, Suburan, Urban: 62dBm/1MHz --> 100MHz
% bs_bw_mhz=100; %%%%%%100 MHz channels for base station
% network_loading_reduction=8.25
% pol_mismatch=0 %%%%%%%%%Polarization mismatch dB
% sim_scale_factor=2; %%%%%%%%%%The scaling Factor for the simulation area.For Macrocell --> 2, For Microcell -->3
% max_itm_dist_km=400; %%%%%It just makes it easy if we have a max number
% mitigation_dB=0:30:30;  %%%%%%%%% in dB%%%%% Beam Muting or PRB Blanking (or any other mitigation mechanism):  30 dB reduction %%%%%%%%%%%%Consider have this be an array, 3dB step size, to get a more granular insight into how each 3dB mitigation reduces the coordination zone.
% reliability=50%
% Tpol=1; %%%polarization for ITM
% FreqMHz=7125;
% confidence=50;
% tx_height_m=25%NaN(1,1); %%%%If NaN, then keep the normal heights: base_station_height
% tf_clutter=0;%1;  %%%%%%This if for P2108.
% p2p_idx=1
% sim_folder1='C:\Local Matlab Data\7GHz Non-Aggregate Ground Sims'  %%%%%%%%%This is where you will create a simulation folder.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Microcell >6m Example
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev=202; %%%%%%P2P Test: Github Example: Microcell 9m Dodgers
% grid_spacing=1;  %%%%km:
% bs_eirp=57%%%%horzcat(57,57,57); %%%%%EIRP [dBm/100MHz] for Rural, Suburan, Urban: 37/1MHz --> 100MHz
% bs_bw_mhz=100; %%%%%%100 MHz channels for base station
% network_loading_reduction=8.3
% pol_mismatch=0 %%%%%%%%%Polarization mismatch dB
% sim_scale_factor=3; %%%%%%%%%%The scaling Factor for the simulation area.For Macrocell --> 2, For Microcell -->3
% max_itm_dist_km=400; %%%%%It just makes it easy if we have a max number
% mitigation_dB=0;  %%%%%%%%% in dB%%%%% Beam Muting or PRB Blanking (or any other mitigation mechanism)  
% reliability=50%
% Tpol=1; %%%polarization for ITM
% FreqMHz=7125;
% confidence=50;
% tx_height_m=9%NaN(1,1); %%%%If NaN, then keep the normal heights: base_station_height
% tf_clutter=0;%1;  %%%%%%This if for P2108.
% p2p_idx=1
% sim_folder1='C:\Local Matlab Data\7GHz Non-Aggregate Ground Sims'  %%%%%%%%%This is where you will create a simulation folder.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Microcell <6m Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev=203; %%%%%%P2P Test: Github Example: Microcell 9m Dodgers
grid_spacing=0.25;  %%%%km: (This is a smaller step size because of additional loss of clutter)
bs_eirp=57%%%%horzcat(57,57,57); %%%%%EIRP [dBm/100MHz] for Rural, Suburan, Urban: 37/1MHz --> 100MHz
bs_bw_mhz=100; %%%%%%100 MHz channels for base station
network_loading_reduction=8.3
pol_mismatch=0 %%%%%%%%%Polarization mismatch dB
sim_scale_factor=3; %%%%%%%%%%The scaling Factor for the simulation area.For Macrocell --> 2, For Microcell -->3
max_itm_dist_km=400; %%%%%It just makes it easy if we have a max number
mitigation_dB=0;  %%%%%%%%% in dB%%%%% Beam Muting or PRB Blanking (or any other mitigation mechanism)  
reliability=50%
Tpol=1; %%%polarization for ITM
FreqMHz=7125;
confidence=50;
tx_height_m=5.9%NaN(1,1); %%%%If NaN, then keep the normal heights: base_station_height
tf_clutter=1;  %%%%%%This if for P2108.
p2p_idx=1
sim_folder1='C:\Local Matlab Data\7GHz Non-Aggregate Ground Sims'  %%%%%%%%%This is where you will create a simulation folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create the cell_sim_data
cell_sim_data=cell_p2p_data(p2p_idx,:);
bs_eirp_reductions=(bs_eirp-bs_down_tilt_reduction-network_loading_reduction); %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Create a Rev Folder
cd(sim_folder1);
pause(0.1)
tempfolder=strcat('Rev',num2str(rev));
[status,msg,msgID]=mkdir(tempfolder);
rev_folder=fullfile(sim_folder1,tempfolder);
cd(rev_folder)
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Saving the simulation files in a folder for the option to run from a server
'First save . . .' 
tic;
save('grid_spacing.mat','grid_spacing')
save('reliability.mat','reliability')
save('confidence.mat','confidence')
save('FreqMHz.mat','FreqMHz')
save('Tpol.mat','Tpol')
save('norm_aas_zero_elevation_data.mat','norm_aas_zero_elevation_data')
save('tf_clutter.mat','tf_clutter')
save('mitigation_dB.mat','mitigation_dB')
save('tx_height_m.mat','tx_height_m')
save('bs_eirp_reductions.mat','bs_eirp_reductions')
toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Should probably pull this out of the loop so we don't have to do it 4,000 times
%%%%%%%%%%Find the ITM Area Pathloss for the distance array
tic;
max_rx_height=max(cell2mat(cell_sim_data(:,5)))
[array_dist_pl]=itm_area_dist_array_sea_rev2(app,reliability,tx_height_m,max_rx_height,max_itm_dist_km,FreqMHz);
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Saving the simulation files in a folder for the option to run from a server
%%%%%%%%%For Loop the Locations
[num_locations,~]=size(cell_sim_data)
base_id_array=1:1:num_locations; %%%%ALL
table([1:num_locations]',cell_sim_data(:,1))

for base_idx=1:1:num_locations
    tic;
    temp_single_cell_sim_data=cell_sim_data(base_idx,:);
    data_label1=temp_single_cell_sim_data{1};
    data_label1=data_label1(find(~isspace(data_label1)));  %%%%%%%%%%Remove the White Spaces
    cell_sim_data{base_idx,1}=data_label1;

    %%%%%%%%%Make a Folder each Location/System
    cd(rev_folder);
    pause(0.1)
    tempfolder2=strcat(data_label1);
    [status,msg,msgID]=mkdir(tempfolder2);
    sim_folder=fullfile(rev_folder,tempfolder2);
    cd(sim_folder)
    pause(0.1)

    %%%%%%%%%%%%%%Calculate FDR --> Required Pathloss --> FSPL base_buffer
    rx_bw_mhz=temp_single_cell_sim_data{4};
    fdr_dB=10*log10(bs_bw_mhz/rx_bw_mhz);
    cell_sim_data{base_idx,13}=fdr_dB;

    rx_nf=temp_single_cell_sim_data{10};
    rx_ant_gain_mb=temp_single_cell_sim_data{9};
    in_ratio=temp_single_cell_sim_data{11};
    min_ant_loss=temp_single_cell_sim_data{12};%     % % 12) Main to side gain: 40dB

    dpa_threshold=-174+10*log10(rx_bw_mhz*10^6)+rx_nf-rx_ant_gain_mb+fdr_dB+in_ratio+pol_mismatch;
    cell_sim_data{base_idx,14}=dpa_threshold;

    required_pathloss=ceil(bs_eirp_reductions-dpa_threshold); %%%%%%%%%%%%%%%%%Round up
    cell_sim_data{base_idx,15}=required_pathloss;

    %%%%%%%%%%%%Might need to check if it's more than 1 set of lat/lon points
    rx_height_m=temp_single_cell_sim_data{5};
    base_protection_pts=horzcat(temp_single_cell_sim_data{:,[2,3]});
    base_protection_pts(:,3)=rx_height_m;
    save(strcat(data_label1,'_base_protection_pts.mat'),'base_protection_pts')

    base_polygon=base_protection_pts;
    save(strcat(data_label1,'_base_polygon.mat'),'base_polygon')

    min_azimuth=temp_single_cell_sim_data{7};
    max_azimuth=temp_single_cell_sim_data{8};
    ant_beamwidth=temp_single_cell_sim_data{6};

    %%%%%%%%Calculate the pathloss as a function of azimuth
    [azi_required_pathloss]=calc_azi_pathloss_rev1(app,ant_beamwidth,min_ant_loss,min_azimuth,max_azimuth,required_pathloss,array_dist_pl,sim_scale_factor);

    %%%%%%%%%%%azi_required_pathloss
    %%%%%%%%%1)Azimuth Degrees, 2) Pathloss, 3) Distance km for base_buffer
    save(strcat(data_label1,'_azi_required_pathloss.mat'),'azi_required_pathloss')

    %%%%%%%%Make the wider_keyhole
    if tf_clutter==1
        clutter_azi_required_pathloss=azi_required_pathloss;
        clutter_azi_required_pathloss(:,2)=clutter_azi_required_pathloss(:,2)-30;
        clutter_nn_idx=nearestpoint_app(app,clutter_azi_required_pathloss(:,2),array_dist_pl(:,2),'next');
        clutter_azi_required_pathloss(:,3)=array_dist_pl(clutter_nn_idx,1)*sim_scale_factor;
        azi_required_pathloss=clutter_azi_required_pathloss;
        under10km_idx=find(azi_required_pathloss(:,3)<10);
        azi_required_pathloss(under10km_idx,3)=10;
    end

    wider_keyhole=azi_required_pathloss;
    wider_min_azimuth=floor(min_azimuth-(5*ant_beamwidth));
    wider_max_azimuth=ceil(max_azimuth+(5*ant_beamwidth));
    if wider_min_azimuth<0
        mod_wider_min_azimuth=mod(wider_min_azimuth,360);
        min1_idx=nearestpoint_app(app,mod_wider_min_azimuth,wider_keyhole(:,1));
        max2_idx=nearestpoint_app(app,wider_max_azimuth,wider_keyhole(:,1));
        max_dist=max(wider_keyhole(:,3));
        wider_keyhole([min1_idx:1:end],3)=max_dist;
        wider_keyhole([1:1:max2_idx],3)=max_dist;
    elseif wider_max_azimuth>360
        mod_wider_max_azimuth=mod(wider_max_azimuth,360);
        min1_idx=nearestpoint_app(app,wider_min_azimuth,wider_keyhole(:,1));
        max2_idx=nearestpoint_app(app,mod_wider_max_azimuth,wider_keyhole(:,1));
        max_dist=max(wider_keyhole(:,3));
        wider_keyhole([min1_idx:1:end],3)=max_dist;
        wider_keyhole([1:1:max2_idx],3)=max_dist;
    else
        min1_idx=nearestpoint_app(app,wider_min_azimuth,wider_keyhole(:,1));
        max2_idx=nearestpoint_app(app,wider_max_azimuth,wider_keyhole(:,1));
        max_dist=max(wider_keyhole(:,3));
        wider_keyhole([min1_idx:1:max2_idx],3)=max_dist;
    end

    save(strcat(data_label1,'_wider_keyhole.mat'),'wider_keyhole')
    strcat(num2str(base_idx/num_locations*100),'%')
    toc;
end

cd(rev_folder)
pause(0.1)
tic;
save('cell_sim_data.mat','cell_sim_data')
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now running the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf_server_status=0;
parallel_flag=0;
tf_recalculate=0;
[workers,parallel_flag]=check_parallel_toolbox(app,parallel_flag);
wrapper_bugsplat_widen_azimuth_rev10(app,rev_folder,parallel_flag,workers,tf_server_status,tf_recalculate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end_clock=clock;
total_clock=end_clock-top_start_clock;
total_seconds=total_clock(6)+total_clock(5)*60+total_clock(4)*3600+total_clock(3)*86400;
total_mins=total_seconds/60;
total_hours=total_mins/60;
if total_hours>1
    strcat('Total Hours:',num2str(total_hours))
elseif total_mins>1
    strcat('Total Minutes:',num2str(total_mins))
else
    strcat('Total Seconds:',num2str(total_seconds))
end
cd(folder1)
'Done'