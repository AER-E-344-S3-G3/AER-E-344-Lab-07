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
L = 0.101; % [m]
% Curve fit equation for voltage (x) vs air speed (y) from Lab 6
v_vs_U = @(x) -4982.7570*x.^4 + 24021.8937*x.^3 + -43230.9413*x.^2 + 34446.6867*x + -10248.1714;
Cf = [34446.6867 -43230.9413 24021.8937 -4982.7570];
% v_vs_U = @(x) (6553.11832)*x^4 + (-30178.43731)*x^3 + (52087.72644)*x^2 + (-39881.78843)*x + (11425.60699);
% Cf = [-39881.78843 52087.72644 -30178.43731 6555.1183];
row = length(aoa);
col = length(pos);
volt = zeros(row, col); % [V] voltage
U = zeros(row, col); % [m/s] velocity
Ue = 1.32372*(15) -0.430966; % [m/s] Freestream velocity estimated from Lab 2
d_m = .2*.0254; % [m] Distance between measurements

%% Reynolds Number
Re = rho_air * Ue * L / mu_air;

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
xticks = get(gca, 'XTick');
abs_xticks = abs(xticks);
set(gca, 'XTickLabel', abs_xticks);
%saveas(gcf, figure_dir + title_str + ".svg");

%% Position vs. Voltage Graph
% figure(figure_count);
% figure_count = figure_count+1;
% for i = 1:row
%     plot(pos, volt(i,:));
%     hold on
% end
% fontname("Times New Roman");
% fontsize(12, "points");
% title_str = "Position vs. Anemometer Voltage at 4° to 16° AOA";
% title(title_str);
% xlabel("Position [in]");
% ylabel("Voltage [V]");
% grid on;
% legend(aoa + "° AOA");

% Boundry Layer determined from graph estimations (looking/making it up)
delta = [0.4, 1, 1.5, 3]; % [in]
delta_pos = 1.8; % [in] try to make position centered at wake
y_right = .2 : .2 : max(pos) - delta_pos;
y_left = min(pos) - delta_pos : .2 : 0;
y_pos = [y_left,y_right];

%% Normalized Position vs. Normalized Velocity Graph
y_delta = zeros(row,col);
U_Ue = zeros(row,col);
for i = 1:row
    figure(figure_count)
    figure_count = figure_count+1;
    y_delta(i,:) = y_pos/delta(i);
    U_Ue(i,:) = U(i,:)/Ue;
    plot(y_delta(i,:), U_Ue(i,:));
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Y/\delta vs. U/Ue at "+ aoa(i)+"° AOA";
    title(title_str);
    xlabel("Y/\delta");
    ylabel("U/Ue");
    grid on;
    xticks = get(gca, 'XTick');
    abs_xticks = abs(xticks);
    set(gca, 'XTickLabel', abs_xticks);
    %saveas(gcf, figure_dir + title_str + ".svg");
end

%% Normalized Position vs. Turbulence Intensity 
turb_i = zeros(row,col);
for i = 1:row
    v_rms = zeros(1,col);
    for j = 1:col
        % find rms of raw data at each position
        name = "./" + aoa(i) + "/" + pos(j) + "in.txt";
        temp = readmatrix(name);
        v = temp(4:1003,2);
        v_rms(j) = rms(v);
    end
    v_ = volt(i,:);
    u_rms = (Cf(1) + 2*Cf(2)*v_ + 3*Cf(3)*v_.^2 + 4*Cf(4)*v_.^3).*v_rms;
    turb_i(i,:) = u_rms./U(i,:);

    figure(figure_count)
    figure_count = figure_count+1;
    plot(y_delta(i,:), turb_i(i,:));
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Y/\delta vs. Turbulence Intensity at "+ aoa(i)+"° AOA";
    title(title_str);
    xlabel("Y/\delta");
    ylabel("Turbulence Intensity");
    grid on;
    xticks = get(gca, 'XTick');
    abs_xticks = abs(xticks);
    set(gca, 'XTickLabel', abs_xticks);
    %saveas(gcf, figure_dir + title_str + ".svg");
end

%% Vortex Shedding Frequency
for i = 1:4
    figure(figure_count)
    figure_count = figure_count+1;
    % Inside Wake
    name = "./" + aoa(i) + "/" + pos(10) + "in.txt";
    temp = readmatrix(name);
    y = temp(4:1003,2);
    Fs = 10000;
    n = length(y);
    Y = fft(y);
    P2 = abs(Y/n); % Compute the two-sided spectrum
    P1 = P2(1:n/2+1); % Compute the single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(n/2))/n; % Define the frequency domain 
    % Plot the single-sided amplitude spectrum inside wake
    plot(f,P1);
    hold on

    % Outside Wake
    name = "./" + aoa(i) + "/" + pos(1) + "in.txt";
    temp = readmatrix(name);
    y = temp(4:1003,2);
    Fs = 10000;
    n = length(y);
    Y = fft(y);
    P2 = abs(Y/n); % Compute the two-sided spectrum
    P1 = P2(1:n/2+1); % Compute the single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(n/2))/n; % Define the frequency domain 
    % Plot the single-sided amplitude spectrum inside wake
    plot(f,P1);

    title_str ="Single-Sided Amplitude Spectrum at "+ aoa(i)+"° AOA";
    title(title_str);
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    ylim([0,.3]);
    legend({'Inside Wake', 'Outisde Wake'})
    %saveas(gcf, figure_dir + title_str + ".svg");
end

%% Momentum Thickness
theta = zeros(1,row);
C_D = zeros(1,row);
for i = 1:row
    for j = 10:20 
        f1 = U(i,j)/Ue*(1- U(i,j)/Ue);
        f2 = U(i,j+1)/Ue*(1-U(i,j+1)/Ue);
        theta(i) = theta(i) + d_m*((f1 + f2)/2);
    end
    C_D(i) = 2 * theta(i) / L;
end

%% Delete unzipped data
dirs = {'4', '8', '12', '16'};
for i = 1:length(dirs)
      rmdir(dirs{i}, 's'); 
end
