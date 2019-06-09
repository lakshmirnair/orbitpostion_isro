function main()

RE = 6371;          % Earth's radius                            [km]
muE = 398600.44;    % Earth gravitational parameter             [km^3/sec^2]
wE = (2*pi/86164);  % Earth rotation velocity aorund z-axis     [rad/sec]
disp(' 3D view  of orbit simulation of satelites');
%values for each corresponding variable
%RAAN=-8.28E-09;
%w=6.30E-01;
%v0=1.64;
%i=9.69E-01;
%a=2.5625E+04;
%e=7.13E-03;

%Above values take from excel sheet 
RAAN=xlsread('new.csv','W2:W2');
w=xlsread('new.csv','AA2:AA2');
v0=xlsread('new.csv','P2:P2');
i=xlsread('new.csv','Y2:Y2');
as=xlsread('new.csv','T2:T2');
e=xlsread('new.csv','R2:R2');
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
step = 1440;  % time step             [s]
t    = t0:step:tf+step;                          % vector of time        [s]

disp(' ');
% DYNAMICS part
cosE0 = (e+cos(v0))./(1+e.*cos(v0));               % cosine of initial eccentric anomaly
sinE0 = (sqrt(1-e^2).*sin(v0))./(1+e.*cos(v0));    %   sine of initial eccentric anomaly
E0 = atan2(sinE0,cosE0);                           % initial eccentric anomaly              [rad]
if (E0<0)                                          % E0?[0,2pi]
    E0=E0+2*pi;
end
tp = (-E0+e.*sin(E0))./n+t0;                       % pass time at the perigee               [s]
M  = n.*(t-tp);                                    % mean anomaly                           [rad]
% Mk = Ek - e*sin(Ek);
% Eccentric anomaly (must be solved iteratively for Ek)
E = zeros(size(t,2),1);
for j=1:size(t,2)
    E(j) = anom_ecc(M(j),e);                     % eccentric anomaly         [rad]
end
% True anomaly, Argument of latitude, Radius
sin_v = (sqrt(1-e.^2).*sin(E))./(1-e.*cos(E));   % sine of true anomaly
cos_v = (cos(E)-e)./(1-e.*cos(E));               % cosine of true anomaly
v     = atan2(sin_v,cos_v);                      % true anomaly              [rad]
theta = v + w;                                   % argument of latitude      [rad]
r     = (a.*(1-e.^2))./(1+e.*cos(v));            % radius                    [km]
% Satellite coordinates
% "Inertial" reference system ECI (Earth Centered Inertial)
xp = r.*cos(theta);                          % In-plane x position (node direction)             [km]
yp = r.*sin(theta);                          % In-plane y position (perpendicular node direct.) [km]
xs = xp.*cos(RAAN)-yp.*cos(i).*sin(RAAN);    % ECI x-coordinate SAT                             [km]
ys = xp.*sin(RAAN)+yp.*cos(i).*cos(RAAN);    % ECI y-coordinate SAT                             [km]
zs = yp.*sin(i);                             % ECI z-coordinate SAT                             [km]
rs = p./(1+e.*cos(theta-w));                 % norm of radius SAT                               [km]
disp("ECI  x coordinates of satelites in the orbit path having time step 1440")
disp(xs);
disp("ECI  y coordinates of satelites in the orbit path having time step 1440")
disp(ys);
disp("ECI  z coordinates of satelites in the orbit path having time step 1440")
disp(zs);



%greenwich0 
%day=xlsread('new.csv','AM3:AM3');
%greenwitch angle calculated by using matlab code which is available in the
%the file titled Untitled.m
greenwich0= 3.592493965498579e+02;
grp = 0;
disp(' 1) View 3D ');  
%disp(' 2) Planisphere (Only Ground Track)');
%disp(' 3) Both');
while ((grp~=1)&&(grp~=2)&&(grp~=3))
    grp =1;
end

if ((grp==1)||(grp==3))
    str = 0;
   
    disp(' View 3D Options ');
     %disp(' 2) Moment of the momentum ');
    disp(' 2) Equator ');
  
    pw  = 2;
    if (exist('pw','var')==0)
        pw = '0';
    end
   
    disp(' 3) Normal mode');
   
    while ((str~=1)&&(str~=2)&&(str~=3))
        str = 3;
    end
end
greenwich0 = greenwich0*pi/180;                
rot_earth  = wE.*(t-t0)+greenwich0;            
for j=1:size(t,2),
    if rot_earth(j) < (-pi)
        nn = ceil(rot_earth(j)/(-2*pi));
        rot_earth(j) = rot_earth(j) + nn*2*pi;
    elseif rot_earth(j) > (pi)
        nn = fix(rot_earth(j)/(2*pi));
        rot_earth(j) = rot_earth(j) - nn*2*pi;
    end
end
LatSSP     = asin(sin(i).*sin(theta));         
LongSSP    = atan2(ys./rs,xs./rs)-rot_earth';   
   
xSSP = RE.*cos(LatSSP).*cos(LongSSP)+1;       
ySSP = RE.*cos(LatSSP).*sin(LongSSP)+1;        
zSSP = RE.*sin(LatSSP)+0.1;                    


image_file = 'land_ocean_ice_2048.jpg';

cdata = imread(image_file);
if ((grp==1)||(grp==3))
    % Plot the equator
    angle_eq = linspace(0,2*pi,361);
    xeq      = (RE*1.0001).*cos(angle_eq);
    yeq      = (RE*1.0001).*sin(angle_eq);
    zeq      = zeros(1,size(angle_eq,2));
    % Open the window for plotting
    figure('Color','k');
    % Definition of axes necessary to rotate the object to whom are reported
    ax = axes('XLim',[-a a],'YLim',[-a a],'ZLim',[-a a],'Color','k');
    % Base point for Earth
    [x,y,z] = ellipsoid(0, 0, 0, RE, RE, RE, 30);
    hold on;
    switch lower(pw)
        case '1'
            view([h1 h2 h3]);
        case '2'
            view(0,0);
        case '3'
            view(120,90);
    end
  
    globe = surface(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
    grid on;
    axis equal;
  
    if ((ra/RE<13.3)&&(str==2))
        if (ra/RE>9)
            zoom(13.3-ra/RE);
        elseif (ra/RE>6.8)
            zoom(10.7-ra/RE);
        elseif (ra/RE>5)
            zoom(9.6-ra/RE);
        elseif (ra/RE>3.5)
            zoom(6.1-ra/RE);
        elseif (ra/RE>2.7)
            zoom(5.4-ra/RE);
        elseif (ra/RE>2.2)
            zoom(5.0-ra/RE);
        elseif (ra/RE>1.5)
            zoom(3.3-ra/RE);
        else
            zoom(2.7-ra/RE);
        end
    elseif (str==3)
        zoom(ra^1.03/RE);
    elseif (str==1)
        zoom(1);
    end
   
    axis([-ra ra -ra ra -ra ra]);
    hg = hgtransform('Parent',ax);
   set(globe,'Parent',hg);
    set(gcf,'Renderer','opengl');
    
    plot3([0,2*ra],[0,0],[0,0],'-.w','LineWidth',1);                  % X = Aries' direction
    plot3([0,0],[0,2*ra],[0,0],'-.w','LineWidth',1);                  % Y
    plot3([0,0],[0,0],[0,2*ra],'-.w','LineWidth',1);                  % Z
    text(2*ra+120,10,0,texlabel('gamma'),'Color','w','FontSize',18);
    text(10,2*ra+120,0,texlabel('Y'),'Color','w','FontSize',10);
    text(0,0,2*ra+140,texlabel('Z'),'Color','w','FontSize',10);    
    plot3(xeq,yeq,zeq,'--w','LineWidth',1);                           % Equator
    plot3([0,2*ra*n1],[0,2*ra*n2],[0,n3],'-.y','LineWidth',1.5);     % Nodes' Line
    text(2*ra*n1-140,2*ra*n2+140,0,texlabel('RAAN'),'Color','y','FontSize',8);    
  
    u = 0;                                  
    for k=1:size(t,2)                       
        if strcmp(pw,'4')                     
            view([xs(k),ys(k),zs(k)]);
        elseif strcmp(pw,'5')
            view([xs(k),ys(k),2]);
        elseif strcmp(pw,'6')
            view(rot_earth(k)*180/pi+89,2);
        elseif strcmp(pw,'0')
            view(120,45);
        end
        Rz = makehgtform('zrotate',rot_earth(k));                          
        X = plot3(xSSP(k),ySSP(k),zSSP(k),'--bo','LineWidth',2,...         
            'MarkerSize',2.5,'MarkerEdgeColor','b','MarkerFaceColor',[0.1,0.8,0.2]);
        set(X,'Parent',hg);           
        set(hg,'Matrix',Rz);         
        drawnow                      
        u = u+1;                      
        XS(k) = plot3(xs(k),ys(k),zs(k),'rv','LineWidth',0.9,...   
                'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor',[0.8,0.2,0.2]);
        pause(0.001);
        if (u>1)
            delete(XS(k-1));               
            XS(k-1) = plot3(xs(k-1),ys(k-1),zs(k-1),'--mx','LineWidth',0.8,'MarkerSize',2); 
            
        end
    end
end

stop = true(size(t,2),1);
if ((grp==2)||(grp==3))
    while (norm(double(stop))~=0)
        for j=1:size(t,2),
            if LongSSP(j) < (-pi)             %
                %nn = fix(-LongSSP(j)/(pi));
                nn = 1;
                LongSSP(j) = LongSSP(j) + 2*nn*pi;
            elseif LongSSP(j) > (pi)
                %nn = fix(LongSSP(j)/(pi));
                nn = 1;
                LongSSP(j) = LongSSP(j) - 2*nn*pi;
            end
            stop(j) = (abs(LongSSP(j))-pi)>0;
        end
    end
    LongSSP = LongSSP.*180/pi;             % Longitude of SSP [deg]
    LatSSP  = LatSSP.*180/pi;              % Latitude  of SSP [deg] 
   
    planisphere = figure('Name','Planisphere');
    hold on;
    set(gca,'XTick',[-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
        15 30 45 60 75 90 105 120 135 150 165 180],'XTickMode','manual');
    set(gca,'YTick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90],'YTickMode','manual');
        
    switch lower(map)
        case '1'
            image_file = 'land_ocean_ice_2048.jpg';
            cdata      = imread(image_file);
            imagesc([-180,180],[90,-90],cdata);
        case '2'
            image_file = 'land_par_mer.jpg';
            cdata      = imread(image_file);
            imagesc([-180,180],[90,-90],cdata);
        otherwise
           
            load('topo.mat','topo','topomap1');               
            load coast;                                       
            plot(long,lat,'b','LineWidth',1.2); 
    end
    grid on;
    axis([-180 180 -90 90]);
    plot([-180 180],[0 0],'--k','LineWidth',1.5);  
    if (trt==1)
        plot(LongSSP(1),LatSSP(1),'go','MarkerSize',2,'MarkerFaceColor',[0 1 0],...  
                    'MarkerEdgeColor','g');
                % cycle for plotting one orbit
        for j=2:size(t,2)  
            if (t(j)-t0<=T+step);
                plot(LongSSP(j),LatSSP(j),'ro','MarkerSize',1,'MarkerFaceColor',[0.8 0.2 0.1],...
                    'MarkerEdgeColor','r');
                pause(0.05);
                if (j>1)
                    if ((abs(LongSSP(j)-LongSSP(j-1))<160)&&(abs(LatSSP(j)-LatSSP(j-1))<135))
                        plot([LongSSP(j-1),LongSSP(j)],[LatSSP(j-1),LatSSP(j)],'-r',...
                            'LineWidth',1.5);
                    end
                end
            end
        end
        plot(LongSSP(end),LatSSP(end),'ko','MarkerSize',2,'MarkerFaceColor',[0 0 0],...   % End
                    'MarkerEdgeColor','w');
    else
        plot(LongSSP(1),LatSSP(1),'go','MarkerSize',2,'MarkerFaceColor',[0 1 0],...   % Departure
                    'MarkerEdgeColor','g');
        for j=2:size(t,2)                 % cycle for plotting all orbits
            plot(LongSSP(j),LatSSP(j),'ro','MarkerSize',1,'MarkerFaceColor',[0.8 0.2 0.2],...
                    'MarkerEdgeColor','r');
            pause(0.05);
            if (j>1)
                if (abs(LongSSP(j)-LongSSP(j-1))<80)
                    plot([LongSSP(j-1),LongSSP(j)],[LatSSP(j-1),LatSSP(j)],'-r',...
                        'LineWidth',1.5);
                end
            end
        end
        plot(LongSSP(end),LatSSP(end),'ko','MarkerSize',2,'MarkerFaceColor',[0 0 0],...   % End
                    'MarkerEdgeColor','w');
    end
end
saveas(gcf,'a1.png');
