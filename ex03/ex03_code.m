
%{Inertial Navigation Lab3 in LUH
% creat by Tingde Liu (1.1.2024)
% E-mail: tingde.liu.luh@gmail.com
%}

%task1 Initial Alignment & Initialisation----------------------------------
load('ex03_data.mat')
w_x = imudata.static(:,2);
w_y = imudata.static(:,3);
w_z = imudata.static(:,4);
f_x = imudata.static(:,5);
f_y = imudata.static(:,6);
f_z = imudata.static(:,7);

%Kresel static (rad/s)
w1 = mean(w_x);
w2 = mean(w_y);
w3 = mean(w_z);

%Beschleunigunsmesser static (m/s^2)
f1 = mean(f_x);
f2 = mean(f_y);
f3 = mean(f_z);

fb_ib = [f1; f2; f3];
wb_ib = [w1; w2; w3];
%initial roll (b-n)
ir = atan2d(-f2, -f3);
display(ir)
%initial pitch (b-n)
ip = atand(f1/ sqrt(f2^2 + f3^2));
display(ip)
%initial yaw (b-n)
C1 = [ 1 0 0 ; 0 cosd(ir) sind(ir) ; 0 -sind(ir) cosd(ir)];
C2 = [ cosd(ip) 0 -sind(ip) ; 0 1 0 ; sind(ip) 0 cosd(ip)];
W_n_ib = C2' * C1' * wb_ib;
iy = atand(-W_n_ib(2)/ W_n_ib(1));
display(iy)

%initial velocity  (n-frame)
an_eb = [0; 0; 0];
vn_eb = [0; 0; 0];
%initial position  (n-frame)
S = [0; 0; 0];

%Earth rotation rate

we = sqrt(W_n_ib(1)^2 + W_n_ib(2)^2 + W_n_ib(3)^2);  % rad/s
we_deg = 3600 * rad2deg(we); % degree/hour
display(we_deg)
%g-vector
g = sqrt(f1^2 + f2^2 + f3^2);  % m/s^2
display(g)


%latitude
wn_be = [0; 0; 0];
wn_ib = W_n_ib;
wn_ie = wn_ib +wn_be;
phi = acos(1 / we * sqrt(wn_ie(1)^2 + wn_ie(2)^2)); % rad
phi_deg = rad2deg(phi); % degree
display(phi_deg)


%Vector
wn_ie =  [we*cos(phi) 0 -we*sin(phi)]';
gn = [0 0 -g]';
% DCM
Cb_n = [fb_ib wb_ib cross(fb_ib, wb_ib)] / [gn wn_ie cross(gn, wn_ie)];
%Reorthogonaliesierung
C_oth = Cb_n * (Cb_n' * Cb_n)^(-0.5);
Cb_n = C_oth;
Cn_b = inv(Cb_n)

%task2 Strapdown Algorithm-------------------------------------------------

%initial value 
an_eb = [0; 0; 0];
vn_eb = [0; 0; 0];
S = [0; 0; 0];
I = [1 0 0;
     0 1 0;
     0 0 1] ;
dt = 0.005; % time interval (data rate:200 Hz)
R0 = 6378000; % earth radius: 6378000m
h = 70; % ellipsoidic height: 70m

Tr = [];    
orientation = [];
velocity = [];
position = [];
beta = [0;0;0];
for i = 1:length(imudata.moving)
    fz = imudata.moving(i,7);
    fy = imudata.moving(i,6);
    fx = imudata.moving(i,5);
    wx = imudata.moving(i,2);
    wy = imudata.moving(i,3);
    wz = imudata.moving(i,4);
    t = imudata.moving(i,1);

    vn_eb_N = vn_eb(1); %Northward velocity in the navigational system
    vn_eb_E = vn_eb(2); %Eastward velocity in the navigational system
    wn_en = [vn_eb_E / (R0+h);             
             -vn_eb_N / (R0+h);
             -vn_eb_E * tan(phi)/ (R0+h)];
    wb_ib = [wx; wy; wz]; % Kreisel
    wb_nb = wb_ib - Cb_n * (wn_en + wn_ie); 

    % Falls Drehraten konstant----------------------------------------
    neu_beta = wb_nb * dt;  % Drehwinkel
    % equivalent Rotation Vector (2st order)
    beta = neu_beta + cross(beta/12,neu_beta);
    beta_norm = norm(beta);
    
    %Skew-symmetric matrix
    beta_schief = [   0      -beta(3)  beta(2)    ; 
                   beta(3)     0       -beta(1)   ;
                   -beta(2)  beta(1)      0     ] ;
    %orientation update (DCM)
    neu_Cn_b = Cn_b * (I + sin(beta_norm)/beta_norm*beta_schief  + (1-cos(beta_norm))/(beta_norm^2)*(beta_schief^2));

    %Reorthogonaliesierung
    neu_Cn_b = neu_Cn_b * (neu_Cn_b' * neu_Cn_b)^(-0.5) ;  
    
    %velocity update
    fb_ib = [fx ; fy ; fz];
    neu_an_eb = Cn_b * fb_ib  - gn - cross((2 * wn_ie + wn_en) , vn_eb); 
    neu_vn_eb = vn_eb + 0.5 * (an_eb + neu_an_eb) * dt;

    %position update
    S = S + 0.5 * (vn_eb + neu_vn_eb) * dt;
    
    an_eb = neu_an_eb;
    vn_eb = neu_vn_eb;
    
    Cn_b = neu_Cn_b;
    Cb_n = inv(Cn_b);
    %neue Euler winkel
    phi_deg = rad2deg(atan2(Cn_b(3,2),Cn_b(3,3)));    
    theta_deg = rad2deg(-asin(Cn_b(3,1)));
    psi_deg = rad2deg(atan2(Cn_b(2,1),Cn_b(1,1)));

    %Result (orientation, velocity, position)
    Euler = [phi_deg; theta_deg; psi_deg;t];
    orientation = [orientation, Euler];
    v = [vn_eb; t];
    s = [ S; t];
    tr = [wn_en; t];
    velocity = [velocity, v];
    position = [position, s];
    Tr = [Tr, tr];
end

%orientierung

figure(1);
x = orientation(4,:);
y1 = orientation(1,:);
y2 = orientation(2,:);
y3 = orientation(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('roll[deg]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('pitch[deg]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('yaw[deg]', fontsize=13)


%velocity

figure(2);
x = velocity(4,:);
y1 = velocity(1,:);
y2 = velocity(2,:);
y3 = velocity(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('velocity_x[m/s]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('velocity_y[m/s]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('velocity_z[m/s]', fontsize=13)

%position

figure(3);
x  = position(4,:);
y1 = position(1,:);
y2 = position(2,:);
y3 = position(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('position_x[m]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('position_y[m]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('position_z[m]', fontsize=13)



x = position(1,:);
y = position(2,:);
figure
plot(x,y)
xlabel('Sx[m]', fontsize=13)
ylabel('Sy[m]', fontsize=13)


%transportation rate

figure
x = Tr(4,:);
y1 = Tr(1,:);
y2 = Tr(2,:);
y3 = Tr(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]')
ylabel('transportation rate_x[1/s]')
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]')
ylabel('transportation rate_y[1/s]')
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]')
ylabel('transportation rate_z[1/s]')


%slope 斜率

for i = 1:length(length(imudata.moving))
    slope(i) = asind(sqrt((sind(orientation(1,i)))^2 + (sind(orientation(2,i)))^2));
end

slope_m = mean(slope);





%Task3 Performance Analysis-----------------------------------------------

%acceleration bias ---------------------------

%g = 9.81;  
fb_ib = [1.5*g/1000; 1.5*g/1000; 1.5*g/1000+g];
dt = 0.005;
abv = [];
abs = [];
v = [1; 0; 0];
s = [0; 0; 0];
for i = 1:length(imudata.moving)
    t = imudata.moving(i,1);
    v_new = (Cb_n*fb_ib+gn)*dt + v;
    s_new = 0.5*(v_new+v)*dt + s;
    v_t = [v_new;  t];
    s_t = [s_new;  t];
    abv = [abv, v_t];
    abs = [abs, s_t];
    v = v_new;
    s = s_new;
end

figure
x  = abv(4,:);
y1 = abv(1,:);
y2 = abv(2,:);
y3 = abv(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('ba-velocity_x[m/s]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('ba-velocity_y[m/s]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('ba-velocity_z[m/s]', fontsize=13)

figure
x  = abs(4,:);
y1 = abs(1,:);
y2 = abs(2,:);
y3 = abs(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('ba-position_x[m]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('ba-position_y[m]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('ba-position_z[m]', fontsize=13)

x = abs(1,:);
y = abs(2,:);
figure
plot(x,y)
subtitle('track with acceleration bias')
xlabel('Sx[m]', fontsize=13)
ylabel('Sy[m]', fontsize=13)

%bg einsetzen--------------------------------------------------------------

fb_ib = [0; 0; g];
wb_ib = [0; 0; (0.75*pi)/(3600*180)];

%wn_ie =  [0 0 0]'; %(neglected Coriolis term)

% DCM

psi_deg = 0;
vn_nb = [1; 0; 0];
vb_nb = [1; 0; 0];
S = [0; 0; 0];
dt = 0.005;
trbv = [];
trbs = [];
for i = 1:length(imudata.moving)
   
    t = imudata.moving(i,1);

    wb_nb = wb_ib; %- Cb_n * (wn_en + wn_ie); 

    beta = wb_nb * dt;      
    psi_deg = beta(3) + psi_deg;
    neu_Cn_b = [cos(psi_deg) sin(psi_deg) 0;-sin(psi_deg) cos(psi_deg) 0;0 0 1];
 

    neu_vn_nb = neu_Cn_b*vb_nb ;
    S = S + 0.5 * (vn_nb + neu_vn_nb) * dt;
    
    vn_nb = neu_vn_nb;
    Cn_b = neu_Cn_b;
    Cb_n = inv(Cn_b);

    vt = [vn_eb; t];
    st = [S; t];
    trbv = [trbv, vt];
    trbs = [trbs, st];
end
figure
x = trbv(4,:);
y1 = trbv(1,:);
y2 = trbv(2,:);
y3 = trbv(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('bg-velocity_x[m/s]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('bg-velocity_y[m/s]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('bg-velocity_z[m/s]', fontsize=13)

figure
x = trbs(4,:);
y1 = trbs(1,:);
y2 = trbs(2,:);
y3 = trbs(3,:);
subplot(3,1,1)
plot(x,y1)
xlabel('time[s]', fontsize=13)
ylabel('bg-position_x[m]', fontsize=13)
subplot(3,1,2)
plot(x,y2)
xlabel('time[s]', fontsize=13)
ylabel('bg-position_y[m]', fontsize=13)
subplot(3,1,3)
plot(x,y3)
xlabel('time[s]', fontsize=13)
ylabel('bg-position_z[m]', fontsize=13)

x = trbs(1,:);
y = trbs(2,:);
figure
plot(x,y)
subtitle('track with turn-rate bias bias', fontsize=13)
xlabel('Sx[m]', fontsize=13)
ylabel('Sy[m]', fontsize=13)

%task4----------------------------------------------------------

w_x = imudata.static(:,2);
w_y = imudata.static(:,3);
w_z = imudata.static(:,4);
f_x = imudata.static(:,5);
f_y = imudata.static(:,6);
f_z = imudata.static(:,7);

%Kresel static (rad/s)
w1 = mean(w_x);
w2 = mean(w_y);
w3 = mean(w_z);

%Beschleunigunsmesser static (m/s^2)
f1 = mean(f_x);
f2 = mean(f_y);
f3 = mean(f_z);

fb_ib = [f1; f2; f3];
wb_ib = [w1; w2; w3];
%initial roll (b-n)
%ir = atan2d(-f2, -f3);
display(ir)
%initial pitch (b-n)
ip = atand(f1/ sqrt(f2^2 + f3^2));
%display(ip)
%initial yaw (b-n)
C1 = [ 1 0 0 ; 0 cosd(ir) sind(ir) ; 0 -sind(ir) cosd(ir)];
C2 = [ cosd(ip) 0 -sind(ip) ; 0 1 0 ; sind(ip) 0 cosd(ip)];
W_n_ib = C2' * C1' * wb_ib;
iy = atand(-W_n_ib(2)/ W_n_ib(1));
%display(iy)

%initial velocity  (n-frame)
an_eb = [0; 0; 0];
vn_eb = [0; 0; 0];
%initial position  (n-frame)
S = [0; 0; 0];

%Earth rotation rate

we = sqrt(W_n_ib(1)^2 + W_n_ib(2)^2 + W_n_ib(3)^2);  % rad/s
we_deg = 3600 * rad2deg(we); % degree/hour
%display(we_deg)
%g-vector
g = sqrt(f1^2 + f2^2 + f3^2);  % m/s^2
%display(g);


%latitude
wn_be = [0; 0; 0];
wn_ib = W_n_ib;
wn_ie = wn_ib +wn_be;
phi = acos(1 / we * sqrt(wn_ie(1)^2 + wn_ie(2)^2)); % rad
phi_deg = rad2deg(phi); % degree
%display(phi_deg);


%Vector
wn_ie =  [we*cos(phi) 0 -we*sin(phi)]';
gn = [0 0 -g]';
% DCM
Cb_n = [fb_ib wb_ib cross(fb_ib, wb_ib)] / [gn wn_ie cross(gn, wn_ie)];
%Reorthogonaliesierung
C_oth = Cb_n * (Cb_n' * Cb_n)^(-0.5);
Cb_n = C_oth;
Cn_b = inv(Cb_n);
%initial value 
an_eb = [0; 0; 0];
vn_eb = [0; 0; 0];
S = [0; 0; 0];
I = [1 0 0;
     0 1 0;
     0 0 1] ;
dt = 0.05; % time interval (data rate:20 Hz)
R0 = 6378000; % earth radius: 6378000m
h = 70; % ellipsoidic height: 70m

Tr = [];    
orientation = [];
velocity = [];
position = [];
beta = [0;0;0];
for i = 1:10:length(imudata.moving)
    fz = imudata.moving(i,7);
    fy = imudata.moving(i,6);
    fx = imudata.moving(i,5);
    wx = imudata.moving(i,2);
    wy = imudata.moving(i,3);
    wz = imudata.moving(i,4);
    t = imudata.moving(i,1);

    vn_eb_N = vn_eb(1); %Northward velocity in the navigational system
    vn_eb_E = vn_eb(2); %Eastward velocity in the navigational system
    wn_en = [vn_eb_E / (R0+h);             
             -vn_eb_N / (R0+h);
             -vn_eb_E * tan(phi)/ (R0+h)];
    wb_ib = [wx; wy; wz]; % Kreisel
    wb_nb = wb_ib - Cb_n * (wn_en + wn_ie); 

    % Falls Drehraten konstant----------------------------------------
    neu_beta = wb_nb * dt;  % Drehwinkel
    % equivalent Rotation Vector (2st order)
    beta = neu_beta + cross(beta/12,neu_beta);
    beta_norm = norm(beta);
    
    %Skew-symmetric matrix
    beta_schief = [   0      -beta(3)  beta(2)    ; 
                   beta(3)     0       -beta(1)   ;
                   -beta(2)  beta(1)      0     ] ;
    %orientation update (DCM)
    neu_Cn_b = Cn_b * (I + sin(beta_norm)/beta_norm*beta_schief  + (1-cos(beta_norm))/(beta_norm^2)*(beta_schief^2));

    %Reorthogonaliesierung
    neu_Cn_b = neu_Cn_b * (neu_Cn_b' * neu_Cn_b)^(-0.5) ;  
    
    %velocity update
    fb_ib = [fx ; fy ; fz];
    neu_an_eb = Cn_b * fb_ib  - gn - cross((2 * wn_ie + wn_en) , vn_eb); 
    neu_vn_eb = vn_eb + 0.5 * (an_eb + neu_an_eb) * dt;

    %position update
    S = S + 0.5 * (vn_eb + neu_vn_eb) * dt;
    
    an_eb = neu_an_eb;
    vn_eb = neu_vn_eb;
    
    Cn_b = neu_Cn_b;
    Cb_n = inv(Cn_b);
    %neue Euler winkel
    phi_deg = rad2deg(atan2(Cn_b(3,2),Cn_b(3,3)));    
    theta_deg = rad2deg(-asin(Cn_b(3,1)));
    psi_deg = rad2deg(atan2(Cn_b(2,1),Cn_b(1,1)));

    %Result (orientation, velocity, position)
    Euler = [phi_deg; theta_deg; psi_deg;t];
    orientation = [orientation, Euler];
    v = [vn_eb; t];
    s = [ S; t];
    tr = [wn_en; t];
    velocity = [velocity, v];
    position = [position, s];
    Tr = [Tr, tr];
end




x = position(1,:);
y = position(2,:);
figure
plot(x,y)
xlabel('Sx[m]', fontsize=13)
ylabel('Sy[m]', fontsize=13)


