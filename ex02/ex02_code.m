load('ex02_task1.mat');
load('ex02_task2.mat');
load('ex02_task3.mat');
load('imubsp_static1h_10hz.mat')

%Task1----------------------------------------------------------------------

%Assumptions:
g =  9.811117;, fontsize=14, fontsize=14, fontsize=14
we = 7.2921151467e-005;
phi = 52.385828;

%accelerometer
fxa = mean(sixpos.fx_up);
fxb = mean(sixpos.fx_down);
fya = mean(sixpos.fy_up);
fyb = mean(sixpos.fy_down);
fza = mean(sixpos.fz_up);
fzb = mean(sixpos.fz_down);

%gyro
wxa = mean(sixpos.omx_up);
wxb = mean(sixpos.omx_down);
wya = mean(sixpos.omy_up);
wyb = mean(sixpos.omy_down);
wza = mean(sixpos.omz_up);
wzb = mean(sixpos.omz_down);

%bias and scale factor of the accelerometer
ba1 = 0.5 * (fxa + fxb)
ba2 = 0.5 * (fya + fyb)
ba3 = 0.5 * (fza + fzb)
sa1 = (fxb - fxa) / (2 * g)
sa2 = (fyb - fya) / (2 * g)
sa3 = (fzb - fza) / (2 * g)

%bias and scale factor of the gyro
bg1 = 0.5 * (wxa + wxb)
bg2 = 0.5 * (wya + wyb)
bg3 = 0.5 * (wza + wzb) 
sg1 = (wxb - wxa) / (2 * we * sind(phi))
sg2 = (wyb - wya) / (2 * we * sind(phi))
sg3 = (wzb - wza) / (2 * we * sind(phi))



%Task2-------------------------------------------------------------------
% a) Plot 

x1 = 1:10000;
x2 = 1:9999;
x3 = 1:9998;

ya = sim_processes.a;
ya1 = diff(sim_processes.a);
ya2 = diff(sim_processes.a,2);

yb = sim_processes.b;
yb1 = diff(sim_processes.b);
yb2 = diff(sim_processes.b,2);

yc = sim_processes.c;
yc1 = diff(sim_processes.c);
yc2 = diff(sim_processes.c,2);

yd = sim_processes.d;
yd1 = diff(sim_processes.d);
yd2 = diff(sim_processes.d,2);

%datei a
%origin

figure(1);
subplot(3,1,1)
set(gca, "FontSize", 28)
plot(x1,ya,'r');
%title('a origin');
xlabel('timeseries', fontsize=14);
ylabel('origin a', fontsize=14);

%differentiation in first order
%figure(2);
subplot(3,1,2)
plot(x2,ya1,'b');
%title('a first diff');
xlabel('timeseries', fontsize=14);
ylabel('diff. in first', fontsize=14);

%differentiation in second order
%figure(3);
subplot(3,1,3)
plot(x3,ya2,'g');
%title('a second diff');
xlabel('timeseries', fontsize=14);
ylabel('diff. in second', fontsize=14);

%datei b
%origin
figure(2);
subplot(3,1,1)
plot(x1,yb,'r');
%title('b origin');
xlabel('timeseries', fontsize=14);
ylabel('origin b', fontsize=14);

%differentiation in first order
subplot(3,1,2)
plot(x2,yb1,'b');
%title('b first');
xlabel('timeseries', fontsize=14);
ylabel('diff. in first', fontsize=14);

%differentiation in second order
subplot(3,1,3)
plot(x3,yb2,'g');
%title('b second');
xlabel('timeseries', fontsize=14);
ylabel('diff. in second', fontsize=14);

%datei c
%origin
figure(3);
subplot(3,1,1)
plot(x1,yc,'r');
%title('c origin');
xlabel('timeseries', fontsize=14);
ylabel('origin c', fontsize=14);

%differentiation in first
%figure(8);
subplot(3,1,2)
plot(x2,yc1,'b');
%title('c first');
xlabel('timeseries', fontsize=14);
ylabel('diff. in first', fontsize=14);

%differentiation in second
%figure(9);
subplot(3,1,3)
plot(x3,yc2,'g');
%title('c second');
xlabel('timeseries', fontsize=14);
ylabel('diff. in second', fontsize=14);

%datei d
%origin
figure(4);
subplot(3,1,1)
plot(x1,yd,'r');
%title('d origin');
xlabel('timeseries', fontsize=14);
ylabel('origin d', fontsize=14);

%differentiation in first
%figure(11);
subplot(3,1,2)
plot(x2,yd1,'b');
%title('d first');
xlabel('timeseries', fontsize=14);
ylabel('diff. in first', fontsize=14);

%differentiation in second
%figure(12);
subplot(3,1,3)
plot(x3,yd2,'g');
%title('d second');
xlabel('timeseries', fontsize=14);
ylabel('diff. in second', fontsize=14);

% b) Compute the Allan-deviation ---------------------------------

%Allan-Deviation

%a
allan_a = [];
for k = 1:fix(10000/6)
ak = sim_processes.a(1,1:(10000 - mod(10000,k)));
a_r = reshape(ak',k,fix(10000/k));
ybar = mean(a_r',2);
ysub = diff(ybar).* diff(ybar);
Sigma_aq = 1 / 2 / (fix(10000/k)-1) * sum(ysub);
Sigma_a = sqrt(Sigma_aq);
allan_a = [allan_a,Sigma_a];
end


%b
allan_b = [];
for k = 1:fix(10000/6)
bk = sim_processes.b(1,1:(10000 - mod(10000,k)));
b_r = reshape(bk',k,fix(10000/k));

ybarb = mean(b_r',2);
ysubb = diff(ybarb).* diff(ybarb);
Sigma_bq = 1 / 2 / (fix(10000/k)-1) * sum(ysubb);
Sigma_b = sqrt(Sigma_bq);
allan_b = [allan_b,Sigma_b];
end


%c
allan_c = [];
for k = 1:fix(10000/6)
ck = sim_processes.c(1,1:(10000 - mod(10000,k)));
c_r = reshape(ck',k,fix(10000/k));

ybarc = mean(c_r',2);
ysubc = diff(ybarc).* diff(ybarc);
Sigma_cq = 1 / 2 / (fix(10000/k)-1) * sum(ysubc);
Sigma_c = sqrt(Sigma_cq);
allan_c = [allan_c,Sigma_c];
end


%d
allan_d = [];
for k = 1:fix(10000/6)
dk = sim_processes.d(1,1:(10000 - mod(10000,k)));
d_r = reshape(dk',k,fix(10000/k));

ybard = mean(d_r',2);
ysubd = diff(ybard).* diff(ybard);
Sigma_dq = 1 / 2 / (fix(10000/k)-1) * sum(ysubd);
Sigma_d = sqrt(Sigma_dq);
allan_d = [allan_d,Sigma_d];
end

%in a commen plot 
x = 1:1666;
y1 = allan_a;
y2 = allan_b;
y3 = allan_c;
y4 = allan_d;
figure(5);
loglog(x,y1,'r');
hold on;
loglog(x,y2,'g');
hold on;
loglog(x,y3,'b');
hold on;
loglog(x,y4,'k');
title('Allan-Deviation', fontsize=14);
xlabel('\tau [s]', fontsize=14);
ylabel('\sigma_y ', fontsize=14);

% c) Estimate the sum-process (a+b+c+d) ----------------------------

abcd = sim_processes.a + sim_processes.b + sim_processes.c + sim_processes.d;
allan_abcd = [];
for k = 1:fix(10000/6)
abcdk = abcd(1,1:(10000 - mod(10000,k)));
abcd_r = reshape(abcdk',k,fix(10000/k));

ybar_abcd = mean(abcd_r',2);
ysub_abcd = diff(ybar_abcd).* diff(ybar_abcd);
Sigmaq_abcd = 1 / 2 / (fix(10000/k)-1) * sum(ysub_abcd);
Sigma_abcd = sqrt(Sigmaq_abcd);
allan_abcd = [allan_abcd,Sigma_abcd];
end
%in a commen plot
x = 1:1666;
y5 = allan_abcd;
figure(6);
loglog(x,y5,'r');
title('Alan-Deviation sum（a,b,c,d）', fontsize=14);
xlabel('\tau[s]', fontsize=14);
ylabel('\sigma_y', fontsize=14);

%Task3-------------------------------------------------------------------\
% a plot the raw data 

t_l = mustrain_10hz(:,1);
wx_l = mustrain_10hz(:,2);
wy_l = mustrain_10hz(:,3);
wz_l = mustrain_10hz(:,4);
fx_l = mustrain_10hz(:,5);
fy_l = mustrain_10hz(:,6);
fz_l = mustrain_10hz(:,7);

x = t_l;
y_wx = wx_l;
y_wy = wy_l;
y_wz = wz_l;
y_ax = fx_l;
y_ay = fy_l;
y_az = fz_l;

figure(7);
subplot(3,1,1)
plot(x,y_wx,'b');
title('The raw data of Gyro(x-axis)', fontsize=14);
xlabel('time[s]', fontsize=14);
ylabel('turn rate[rad/s]', fontsize=14);
subplot(3,1,2)
plot(x,y_wy,'r');
title('The raw data of Gyro(y-axis)', fontsize=14);
xlabel('time[s]', fontsize=14);
ylabel('turn rate[rad/s]', fontsize=14);
subplot(3,1,3)
plot(x,y_wz,'g');
title('The raw data of Gyro(z-axis)', fontsize=14);
xlabel('time[s]', fontsize=14);
ylabel('turn rate[rad/s]', fontsize=14);

figure(8);
subplot(3,1,1)
plot(x,y_ax,'b');
title('The raw data of accelerator(x-axis)');
xlabel('time[s]', fontsize=14);
ylabel('accelerator[m/s**2]', fontsize=14);
subplot(3,1,2)
plot(x,y_ay,'r');
title('The raw data of accelerator(y-axis)');
xlabel('time[s]', fontsize=14);
ylabel('accelerator[m/s**2]', fontsize=14);
subplot(3,1,3)
plot(x,y_az,'g');
title('The raw data of accelerator(z-axis)');
xlabel('time[s]', fontsize=14);
ylabel('accelerator[m/s**2]', fontsize=14);

%IMU zusammen rechnen 

x = t_l;
y_1 = wx_l;
y_2 = wy_l;
y_3 = wz_l;
y_4 = fx_l;
y_5 = fy_l;
y_6 = fz_l;

A = [y_1 y_2 y_3 y_4 y_5 y_6];  
for j = 1:6
    A_j = A(:,j);
    al_j = [];
    for k = 1:fix(36000/6)
    yj_lk = A_j(1:(36000 - mod(36000,k)),1);
    yj_r = reshape(yj_lk,k,fix(36000/k));
    Ysubwj = diff(mean(yj_r',2)).* diff(mean(yj_r',2));
    Sigmaq_yj = 1 / 2 / (fix(36000/k)-1) * sum(Ysubwj);
    Sigma_yj = sqrt(Sigmaq_yj);
    al_j = [al_j,Sigma_yj];
    end
    Datei{j} = al_j;
end


for j = 1:3
    x = 0.1:0.1:600;
    y = Datei{j};
    figure;
    loglog(x,y,'b');
    xlabel('tau [s]', fontsize=14);
    ylabel('sigma_y [rad/s]', fontsize=14);
end

for j = 4:6
    x = 0.1:0.1:600;
    y = Datei{j};
    figure;
    loglog(x,y,'b');
    xlabel('tau [s]', fontsize=14);
    ylabel('sigma_y [m/s^2]', fontsize=14);
end

%Task 4 ------------------------------------------------------------

% find h_2, h_1, h_0  in acc

tau_2a = 400;
Sigma_2a = 0.000366188;
h_2a = 6/(4 * 3.14159^2) * Sigma_2a^2 / tau_2a
Sigma_1a = 0.000139173;
h_1a = Sigma_1a^2 / 2 * log(2)
tau_0a = 1;
Sigma_0a = 0.000384591;
h_0a = 2 * tau_0a * Sigma_0a^2
S_ya = h_2a * 10^(-2) + h_1a * 10^(-1) + h_0a * 10^0;



% find h_2, h_1, h_0 in gyro

tau_2g = 400;
Sigma_2g = 1.51465e-5;
h_2g = 6/(4 * 3.14159^2) * Sigma_2g^2 / tau_2g
Sigma_1g = 0.000010559;
h_1g = Sigma_1g^2 / 2 * log(2)
tau_0g = 1;
Sigma_0g = 0.0000279921;
h_0g = 2 * tau_0g * Sigma_0g^2
S_yg = h_2g * 10^(-2) + h_1g * 10^(-1) + h_0g * 10^0


% FFT 

x = mustrain_10hz(:,4);
fs = 10;
[freq,psdx] = aux_calcPSD(x, fs);
x_Fg = freq;
y_Fg = psdx;

% FFT :acc

x = mustrain_10hz(:,7);
fs = 10;
[freq,psdx] = aux_calcPSD(x, fs);
x_fa = freq;
y_fa = psdx;
f = freq;
sya_f = h_2a * f.^(-2) + h_1a * f.^(-1) + h_0a * f.^0;
x_via_a = f;
y_via_a = sya_f;
figure(15);
loglog(x_fa,y_fa,'b',x_via_a,y_via_a,'r');
title('PSD bei der Beschleunigung', fontsize=14);
xlabel('f(Hz)', fontsize=16);
ylabel('|P_acc(f)|', fontsize=16);

% FFT : gyro
syg_f = h_2g * f.^(-2) + h_1g * f.^(-1) + h_0g * f.^0;
x_via_g = f;
y_via_g = syg_f;
figure(16);
loglog(x_Fg,y_Fg,'b',x_via_g,y_via_g,'r');
title('PSD bei dem Gyroscope', fontsize=14);
xlabel('f(Hz)', fontsize=16);
ylabel('|P_gyro(f)|', fontsize=16);

