load('ex01_10050301.mat')
%Task4--------------------------------------------------------------------------------------------------

readdate=importdata('ex01_10050301.mat'); 
averagedate=mean(readdate);      % der Durchschnitt jeder Spalte
w_x=averagedate(2);              % der Durchschnitt der Erdrotationsrate in x-Achse
w_y=averagedate(3);              % der Durchschnitt der Erdrotationsrate in y-Achse
w_z=averagedate(4);              % der Durchschnitt der Erdrotationsrate in z-Achse
f_x=averagedate(5);              % der Durchschnitt der spezifischen Kräfte in x-Achse
f_y=averagedate(6);              % der Durchschnitt der spezifischen Kräfte in y-Achse
f_z=averagedate(7);              % der Durchschnitt der spezifischen Kräfte in z-Achse

%Method (Direct calculation of the DCM)
f_ib_b=[f_x;f_y;f_z];             
w_ib_b=[w_x;w_y;w_z]; 
g=norm([f_x,f_y,f_z]);
phi=deg2rad(52.385828);
w_e=norm([w_x,w_y,w_z]);
g__n=[0;0;-g];

w_ie_n=[w_e*cos(phi);0;-w_e*sin(phi)];

c_n_b=[f_ib_b,w_ib_b,cross(f_ib_b,w_ib_b)]/[g__n,w_ie_n,cross(g__n,w_ie_n)];
if c_n_b*c_n_b'~=1
    disp("Direct calculation of the DCM ist not orthogonal.")
end
[U,S,V]=svd(c_n_b);
r_c_n_b=U*V';
r_c_b_n=r_c_n_b';
c_n_b=r_c_n_b;
c_b_n=r_c_b_n;
disp(['c_n_b with method Direct calculation of the DCM:'])
disp(c_n_b)

%Method (calculated from Eula angle)
ir=atan2d(-f_y,-f_z);            % initial roll
ip=atand(f_x/sqrt(f_y^2+f_z^2)); % initial pitch
iy=atand(-w_y/w_x);              % initial yaw
c_b_n=[cosd(ip)*cosd(iy)  -cosd(ir)*sind(iy)+sind(ir)*sind(ip)*cosd(iy) sind(ir)*sind(iy)+cosd(ir)*sind(ip)*cosd(iy);
       cosd(ip)*sind(iy)  cosd(ir)*cosd(iy)+sind(ir)*sind(ip)*sind(iy)  -sind(ir)*cosd(iy)+cosd(ir)*sind(ip)*sind(iy)
       -sind(ip)         sind(ir)*cosd(ip)                          cosd(ir)*cosd(ip)]; % die Rotationsmatrix
c_n_b=c_b_n'; % Matrix-Inversion
disp(['c_n_b with method calculated from Eula angle:'])
disp(c_n_b)
phi=atan(-w_z/w_x);
w_ie_n=[w_e*cos(phi);0;-w_e*sin(phi)];
%Task5------------------------------------------------------------------------------------------------
roll=[];
pitch=[];
yaw=[];

for index = 1:60001                
     w_x=readdate(index,2);                    
     w_y=readdate(index,3);              
     w_z=readdate(index,4);              
     f_x=readdate(index,5);              
     f_y=readdate(index,6);              
     f_z=readdate(index,7);
      

     w_ib_b=[w_x;w_y;w_z]; 
     w_en_n=[0;0;0];
     w_nb_b=w_ib_b-c_n_b*(w_en_n+w_ie_n);

     delta=w_nb_b*0.01;
     n_delta=norm(delta);
     o_delta=[0 -delta(3) delta(2);
              delta(3) 0 -delta(1);
              -delta(2) delta(1) 0];

     r_c_b_n=c_b_n*(eye(3,3)+(sin(n_delta)/n_delta)*o_delta+((1-cos(n_delta))/n_delta^2)*o_delta^2);
     r_c_n_b=r_c_b_n';
     %[U,S,V]=svd(new_c_n_b);
     %r_c_n_b=U*V';
     %r_c_b_n=r_c_n_b';

     a_phi=rad2deg(atan(r_c_b_n(3,2)/r_c_b_n(3,3)));
     a_theta=rad2deg(-atan2(r_c_b_n(3,1),sqrt(1-r_c_b_n(3,1)^2)));
     a_thera_d=-asin(r_c_b_n(3,1));
     a_psi=rad2deg(acos(r_c_b_n(1,1)/cos(a_thera_d)));
     %a_psi=rad2deg(-atan(r_c_b_n(2,1)/r_c_b_n(1,1)));

     c_b_n=r_c_b_n;
     c_n_b=r_c_n_b;


     roll=[roll a_phi];
     pitch=[pitch a_theta];
     yaw=[yaw a_psi];
end

time=readdate(:,1);
column_roll=roll.';
column_pitch=pitch.';
column_yaw=yaw.';
figure(1)
plot(time,column_roll)
figure(2)
plot(time,column_pitch)
figure(3)
plot(time,column_yaw)

