% AER E 344 Spring 2024 Lab 04 Analysis
% Section 3 Group 3
clear, clc, close all;

figure_dir = "../Figures/";
data_zip = "./data.zip";

%% Parameters
aoa = 4:4:16;
pos = 0:0.2:4;
figure_count = 1;

%% Constants
% Viscosity calculated at 21.2°C
% https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601
% .html
mu_air = 18.18e-6; % [Pa*s]
rho_air = 1.225; % [kg/m^3]
C = 0.101; % [m]
% Curve fit equation for voltage (x) vs air speed (y) from Lab 6
v_vs_U = @(x) -4982.7570*x.^4 + 24021.8937*x.^3 + -43230.9413*x.^2 + 34446.6867*x + -10248.1714;
Cf = [34446.6867 -43230.9413 24021.8937 -4982.7570];
row = length(aoa);
col = length(pos);
volt = zeros(row, col); % [V] voltage
U = zeros(row, col); % [m/s] velocity
Ue = 1.32372*(15) -0.430966; % Freestream velocity estimated from Lab 2

%% Find avg. Voltage and Air Speed
unzip(data_zip);
for i = 1 : row
    for j = 1 : col
        name = "./" + aoa(i) + "/" + pos(j) + "in.txt";
        temp = readmatrix(name);
        volt(i, j) = mean(temp(4:1003,2));
        U(i, j) = v_vs_U(volt(i, j));
    end
end 

%% Position vs. Velocity Graph
figure(figure_count)
figure_count = figure_count+1;
for i = 1:row
    plot(pos, U(i,:));
    hold on
end
fontname("Times New Roman");
fontsize(12, "points");
title_str = "Position vs. Air Velocity at 4° to 16° AOA";
title(title_str);
xlabel("Position [in]");
ylabel("Velocity [m/s]");
grid on;
legend(aoa + "° AOA");

%% Position vs. Voltage Graph
figure(figure_count);
figure_count = figure_count+1;
for i = 1:row
    plot(pos, volt(i,:));
    hold on
end
fontname("Times New Roman");
fontsize(12, "points");
title_str = "Position vs. Anemometer Voltage at 4° to 16° AOA";
title(title_str);
xlabel("Position [in]");
ylabel("Voltage [V]");
grid on;
legend(aoa + "° AOA");

% Boundry Layer determined from graph estimations (looking/making it up)
delta = [0.4, 1, 1.5, 3]; % [in]
delta_pos = 1.8; % [in] try to make position centered at wake

%% Normalized Position vs. Normalized Velocity Graph
y_delta = zeros(row,col);
U_Ue = zeros(row,col);
for i = 1:row
    figure(figure_count)
    figure_count = figure_count+1;
    y_delta(i,:) = pos/delta(i);
    U_Ue(i,:) = U(i,:)/Ue;
    plot(y_delta(i,:), U_Ue(i,:));
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Y/\delta vs. U/Ue at "+ aoa(i)+"° AOA";
    title(title_str);
    xlabel("Y/\delta");
    ylabel("U/Ue");
    grid on;
end

%% Normalized Position vs. Turbulence Intensity 
for i = 1:row
    v_rms = rms(volt(i,:));
    v_ = volt(i,:);
    u_ = v_vs_U(v_);
    u_rms = (Cf(1) + 2*Cf(2)*v_ + 3*Cf(3)*v_.^2 + 4*Cf(4)*v_.^3)*v_rms;
    turb_i = u_rms./u_;

    figure(figure_count)
    figure_count = figure_count+1;
    plot(y_delta(i,:), turb_i);
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Y/\delta vs. Turbulence Intensity at "+ aoa(i)+"° AOA";
    title(title_str);
    xlabel("Y/\delta");
    ylabel("Turbulence Intensity");
    grid on;
end

%% Vortex Shedding Frequency

dirs = {'4', '8', '12', '16'};  % Replace with your directories
for i = 1:length(dirs)
      rmdir(dirs{i}, 's'); 
end
