% ==================================================================================
% This demo aims to simulate the high-frequency dynamics of bottom DO using a 1D ecosystem model
% It is also part of Wu et al. (2024, High-frequency dynamics of bottom dissolved oxygen in temperate 
% shelf seas: the joint role of tidal mixing and sediment oxygen demand, Limnol. Oceanogr. Under Review)  
%
% Author Info: Wenfan Wu, PhD, Virginia Institute of Marine Science, William & Mary | wwu@vims.edu
% ==================================================================================
%% Step-1: Model Config
clc;clearvars
Mobj.nml_file = 'D:\WorkDisk\Code-repository\GitHub_Projects\1D_Ecosystem_Model_for_LO_paper\cosine.nml';   % namelist file for model configuration
Mobj.latitude = 38.78;   % Latitute of the simulation site.

Mobj.nyears = 1;            % simulation period (yr)
Mobj.dt =  15/24/60;       % time interval (d) 
Mobj.dz = 0.1;                % depth interval (m)

Mobj.depmax = 20;  % max depth (m)
Mobj.nLevs = Mobj.depmax/Mobj.dz+1;   % # of vertical layers
Mobj.depth = (0:Mobj.dz:Mobj.depmax)';   % vertical layers

% Index fo dissolved oxygen in the model.
iDOX = 11;  % oxygen
                                                     
nTimes = round(365*Mobj.nyears/Mobj.dt);         % # of time steps
time_steps = (1:nTimes)*Mobj.dt;

%% Step2: Load data (Add the model code into your path first)
load('InitCnd_demo.mat')  % initial conditions and atmospheric forcing
load('temp_prof_from_schism.mat');  % temperature profile from the hydrological model (SCHISM).
load('cmap_data.mat')   % colormaps

%% Step-3: Model Run
[TA, DA, diffusivity] = CoSiNE_1D_main(Mobj, InitCnd);

%% Step-4: Time-depth plot of diffusivity
figure('Units', 'normalized', 'Position', [0.22, 0.34, 0.47, 0.24], 'Color', 'w')
pcolor(time_steps, -Mobj.depth, diffusivity/60/60/24)
shading flat
colorbar 
colormap(cmap_lansey)
caxis([0 2]*1e-5) %#ok<*CAXIS>
xlabel('Days')
ylabel('Depth (m)')
title('Diffusivity (m^2/s)')
set(gca, 'Layer', 'Top')

%% Step-5: Visualize the model results (the same as Figure 7 in Wu et al. 2024)
% Time series of bottom DO and its high-frequency signal
time = InitCnd.time;
ind_lev = Mobj.nLevs;
DO_bot = squeeze(TA(ind_lev, iDOX, :))*0.032;
DO_hf = DO_bot-movmean(DO_bot, fix(1/Mobj.dt));

figure('Units', 'normalized', 'Position', [2.021484e-01,2.984375e-01,2.496745e-01,5.281250e-01], 'Color', 'w')
tiledlayout(6,1,'TileSpacing','tight')
nexttile([2 1])
plot(time, DO_bot, 'LineWidth', 1, 'Color', 'k')
box on; grid on
hold on
plot(time, 2*ones(size(time)), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2)
xlim([time(1) time(end)])
ylim([0 16])
yticks(0:4:16)
set(gca, 'FontSize', 10, 'LineWidth', 1)
xticks(time(1):calmonths(1):time(end))
xticklabels(datestr(time(1):calmonths(1):time(end), 'm')) %#ok<*DATST>
ylabel('Bottom DO (mg/L)', 'FontWeight','bold', 'FontSize', 11)
text(323.1466668393783,3.374233128834356,1.4210854715202e-14, '2 mg/L', 'Color', 'k')

nexttile([2 1])
colororder([0.0118 0.2627 0.8745; 0.9765 0.4510 0.0235])
yyaxis left
plot(time, DO_hf, 'LineWidth', 0.25, 'Color', [0.0118 0.2627 0.8745])
ylabel('HF DO signal (mg/L)', 'FontWeight','bold')
ylim([-0.5 0.5])
yticks(-0.5:0.2:0.5)
yyaxis right
plot(Hydro.time, mean(Hydro.temp_prof(1,:), 1, 'omitnan')-mean(Hydro.temp_prof(21,:), 1, 'omitnan'), 'LineWidth', 2, 'Color', [0.9765 0.4510 0.0235])
box on; grid on
set(gca, 'FontSize', 10, 'LineWidth', 1)
ylabel('Thermal Stratification (℃)', 'FontWeight','bold', 'FontSize', 11)
ylim([-15 15])
yticks(-15:5:15)
xlim([time(1) time(end)])
xticks(time(1):calmonths(1):time(end))
xticklabels(datestr(time(1):calmonths(1):time(end), 'm')) 

nexttile([2 1])
colororder([0.0235 0.7608 0.6745; 0.8118    0.3843    0.4588])
yyaxis left
dox_data = squeeze(TA(:, iDOX, :));
dox_grad = -diffxy(Mobj.depth, dox_data, 1);
plot(time, mean(dox_grad(end-10:end, :), 1), 'LineWidth', 0.25, 'Color', [0.0235 0.7608 0.6745])
xlim([time(1) time(end)])
ylim([-5 30])
yticks(-5:5:30)
box on; grid on
set(gca, 'FontSize', 10, 'LineWidth', 1)
xticks(time(1):calmonths(1):time(end))
xticklabels(datestr(time(1):calmonths(1):time(end), 'm')) %#ok<*DATST>
ylabel('Vertical DO gradient (mg/L/m)', 'FontWeight','bold', 'FontSize', 11)
yyaxis right
plot(time, -DA.o2flx(1,:), 'Linewidth', 2, 'Color', [0.8118 0.3843 0.4588])
ylabel('SOD rate (mmol/m^2/d)', 'FontWeight','bold')
ylim([-5 20])

figure('Units', 'normalized', 'Position',  [0.216,0.09,0.38,0.14], 'Color', 'w')
area(time, DO_hf, 'FaceColor', [0.0039 0.3961 0.9882])
xlim([datetime(2021,7,1) datetime(2021,7,10)])
ylim([-0.4 0.4])
yticks(-0.4:0.2:0.4)
set(gca, 'FontSize', 10, 'LineWidth', 1)
xticks(datetime(2021,7,1:10))
xticklabels(datestr(datetime(2021,7,1:10), 'mmm/dd')) %#ok<*DATST>
ylabel('HF DO signal (mg/L)', 'FontWeight','bold', 'FontSize', 11)

% Vertical profile of DO on July 3 as an example
cyear = year(InitCnd.time(1));
[~, ind_time] = min(abs(time-datetime(cyear,7,3)));
DO_bot = squeeze(TA(:, iDOX, ind_time))*0.032;

figure('Units', 'normalized', 'Position', [4.544271e-01,4.770833e-01,1.884766e-01,4.156250e-01], 'Color', 'w')
plot(DO_bot, -Mobj.depth, 'LineWidth', 3, 'Color', [0.0039 0.4824 0.5725])
hold on
xlim([2 16])
ylim([-Mobj.depmax 0])
box on; grid on
set(gca, 'FontSize', 12, 'LineWidth', 1)
xlabel('DO (mg/L)', 'FontWeight','bold', 'FontSize', 14)
ylabel('Depth (m)', 'FontWeight','bold', 'FontSize', 14)

ax2 = axes('Position',[0.35,0.45,0.35,0.4]);
plot(DO_bot, -Mobj.depth, 'LineWidth', 1, 'Color', [0.1059 0.1412 0.1922], 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
set(ax2, 'Color', [.9 .9 .9])
hold on
ylim([-Mobj.depmax 2-Mobj.depmax])
box on; grid on
set(gca, 'FontSize', 10)
xlabel('DO (mg/L)', 'FontWeight','bold', 'FontSize', 12)
ylabel('Depth (m)', 'FontWeight','bold', 'FontSize', 12)

% Scatter-density diagram of high-frequency DO signals against tidal diffusivity 
time_days = (1:nTimes)*Mobj.dt;

freq_12_hours = 1/12.42;
signal_12_hours = 1*sin(2*pi*freq_12_hours*time_days*24);
freq_8_hours = 1/8.32;
signal_8_hours = 1*sin(2*pi*freq_8_hours*time_days*24);
combined_signal = signal_8_hours+signal_12_hours;
diff_bot = 0.00001*normalize(combined_signal, 'range', [0.1, 10]);

timespan = [datetime(2021,5,1) datetime(2021,9, 1)];  % stratified period
[~, ind_beg] = min(abs(time-timespan(1)));
[~, ind_end] = min(abs(time-timespan(2)));

dx = 0.25e-5;
dy =  6000*dx;

diff_bot_stra = diff_bot(ind_beg:ind_end);
do_hf_stra = DO_hf(ind_beg:ind_end);
[counts,xbins,ybins] = scatter_counts(diff_bot_stra(:), do_hf_stra(:), [dx dy]);

figure('Units', 'normalized', 'Position', [4.322917e-01,2.463542e-01,2.281901e-01,2.817708e-01], 'Color', 'w')
pcolor(xbins, ybins, counts)
shading flat
hold on
xlim([0 10e-5])
ylim([-0.2 0.2])
yticks(-0.2:0.1:0.2)
hold on
set(gca, 'Layer', 'top')
set(gca, 'FontSize', 12, 'LineWidth', 1)
colorbar
colormap(cmap_haxby)
caxis([0 90])
box on; grid on
ylabel('Bottom DO (mg/L)', 'FontWeight','bold', 'FontSize', 14)
xlabel('Tidal diffusivity (m^2/s)', 'FontWeight','bold', 'FontSize', 14)
text(5.9e-5, -0.17, 0, ['N = ', num2str(numel(do_hf_stra))], 'FontSize', 18)
%% END













