function main()

RE = 6371;          % Earth's radius                            [km]
muE = 398600.44;    % Earth gravitational parameter             [km^3/sec^2]
wE = (2*pi/86164);  % Earth rotation velocity aorund z-axis     [rad/sec]
disp(' 3D view  of orbit simulation of satelites');
%RAAN=-8.28E-09;
%w=6.30E-01;
%v0=1.64;
%i=9.69E-01;
%a=2.5625E+04;
%e=7.13E-03;

RAAN=xlsread('new.csv','AB3:AB3');
w=xlsread('new.csv','AA3:AA3');
v0=xlsread('new.csv','P3:P3');
i=xlsread('new.csv','Y3:Y3');
as=xlsread('new.csv','T3:T3');
e=xlsread('new.csv','R3:R3');
a=(as*as)/1000;

rp = a*(1-e);               % radius of perigee             [km]
ra = a*(1+e);               % radius of apogee              [km]
Vp = sqrt(muE*(2/rp-1/a));  % velocity at the perigee       [km/s]
Va = sqrt(muE*(2/ra-1/a));  % velocity at the  apogee       [km/s]
n  = sqrt(muE./a^3);        % mean motion                   [rad/s]1
p  = a*(1-e^2);             % semi-latus rectus             [km]
T  = 2*pi/n;                % period                        [s]
h  = sqrt(p*muE);           % moment of the momentum        [km^2/s]
h1 = sin(i)*sin(RAAN);      % x-component of unit vector h
h2 = -sin(i)*cos(RAAN);     % y-component of unit vector h
h3 = cos(i);                % z-component of unit vector h
n1 = -h2/(sqrt(h1^2+h2^2)); % x-component of nodes' line
n2 =  h1/(sqrt(h1^2+h2^2)); % y-component of nodes' line
n3 = 0;                     % z-component of nodes' line
N  = [n1,n2,n3];            % nodes' line (unit vector)

hours   = floor(T/3600);                   % hours   of the orbital period
minutes = floor((T-hours*3600)/60);        % minutes of the orbital period
seconds = floor(T-hours*3600-minutes*60);  % seconds of the orbital period
fprintf('\n Radius of perigee [%10.3f km]       Altitude of perigee [%10.3f km]',rp,rp-RE);
fprintf('\n Radius of  apogee [%10.3f km]       Altitude of  apogee [%10.3f km]',ra,ra-RE);
fprintf('\n Velocity at the perigee [%6.4f km/s]   Velocity at the apogee [%6.4f km/s]',Vp,Va);
fprintf('\n Orbital Period    [%3d h: %3d m: %3d s] ',hours,minutes,seconds);
fprintf('   = [%10.2f s]\n',T);

norb = 1
t0   = 0;                                        % initial time          [s]
tf   = norb*T;                                   % final   time          [s]  
step = 7200  % time step             [s]
t    = t0:step:tf+step;   
disp(t)
disp(zeros(size(t,2),1));