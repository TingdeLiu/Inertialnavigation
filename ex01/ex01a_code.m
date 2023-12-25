load('ex01_10050301.mat');
%Task1
readdate=importdata("ex01_10050301.mat"); 
averagedate=mean(readdate);      % der Durchschnitt jeder Spalte
w_x=averagedate(2);              % der Durchschnitt der Erdrotationsrate in x-Achse
w_y=averagedate(3);              % der Durchschnitt der Erdrotationsrate in y-Achse
w_z=averagedate(4);              % der Durchschnitt der Erdrotationsrate in z-Achse
f_x=averagedate(5);              % der Durchschnitt der spezifischen Kräfte in x-Achse
f_y=averagedate(6);              % der Durchschnitt der spezifischen Kräfte in y-Achse
f_z=averagedate(7);              % der Durchschnitt der spezifischen Kräfte in z-Achse
ir=atan2d(-f_y,-f_z);            % initial roll
ip=atand(f_x/sqrt(f_y^2+f_z^2)); % initial pitch
iy=atand(-w_y/w_x);              % initial yaw
C_b_n=[cosd(ip)*cosd(iy)  -cosd(ir)*sind(iy)+sind(ir)*sind(ip)*cosd(iy) sind(ir)*sind(iy)+cosd(ir)*sind(ip)*cosd(iy);
       cosd(ip)*sind(iy)  cosd(ir)*cosd(iy)+sind(ir)*sind(ip)*sind(iy)  -sind(ir)*cosd(iy)+cosd(ir)*sind(ip)*sind(iy)
       -sind(ip)         sind(ir)*cosd(ip)                          cosd(ir)*cosd(ip)]; % die Rotationsmatrix
C_n_b=C_b_n'; % Matrix-Inversion
disp(['ir:',num2str(ir)])
disp(['ip:',num2str(ip)])
disp(['iy:',num2str(iy)])
disp(['C_n_b:'])
disp(C_n_b)
%Task2
format long
g=sqrt(f_x^2+f_y^2+f_z^2);                  % die lokale Gravitation
w_e=sqrt(w_x^2+w_y^2+w_z^2);  % die Erdrotation sqrt(w_x^2+w_y^2+w_z^2)*(180/pi)*3600
disp(['g:',num2str(g)])
disp(['w_e:',num2str(w_e)])
%Task3
C_b_n=[cosd(ip)*cosd(iy)  -cosd(ir)*sind(iy)+sind(ir)*sind(ip)*cosd(iy) sind(ir)*sind(iy)+cosd(ir)*sind(ip)*cosd(iy);
       cosd(ip)*sind(iy)  cosd(ir)*cosd(iy)+sind(ir)*sind(ip)*sind(iy)  -sind(ir)*cosd(iy)+cosd(ir)*sind(ip)*sind(iy)
       -sind(ip)         sind(ir)*cosd(ip)                          cosd(ir)*cosd(ip)]; % die Rotationsmatrix
C_n_b=C_b_n'; % Matrix-Inversion
de=[C_n_b(3,2)-C_n_b(2,3);C_n_b(1,3)-C_n_b(3,1);C_n_b(2,1)-C_n_b(1,2)];
u_nb=de/norm(de);   % das Drehvektor u
delta=acos((C_n_b(1,1)+C_n_b(2,2)+C_n_b(3,3)-1)/2); % der Rotationswinkel um u
a=sqrt(1+C_n_b(1,1)+C_n_b(2,2)+C_n_b(3,3))/2;
q=[a;(C_n_b(3,2)-C_n_b(2,3))/(4*a);(C_n_b(1,3)-C_n_b(3,1))/(4*a);(C_n_b(2,1)-C_n_b(1,2))/(4*a)]; % das Quaternion
disp(['u_nb:'])
disp(u_nb)
disp(['delta:',num2str(delta)])
disp(['q'])
disp(q)

