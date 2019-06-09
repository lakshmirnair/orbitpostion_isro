%Find the Julian Date of the previous midnight, JD0
JD0 = 3;
JD=17284;
JDmin = floor(JD)-.5;
JDmax = floor(JD)+.5;
JD0(JD > JDmin) = JDmin(JD > JDmin);
JD0(JD > JDmax) = JDmax(JD > JDmax);
H = (JD-JD0).*24;       %Time in hours past previous midnight
D = JD - 2451545.0;     %Compute the number of days since J2000
D0 = JD0 - 2451545.0;   %Compute the number of days since J2000
T = D./36525;           %Compute the number of centuries since J2000
%Calculate GMST in hours (0h to 24h) ... then convert to degrees
GMST = mod(6.697374558 + 0.06570982441908.*D0  + 1.00273790935.*H + ...
    0.000026.*(T.^2),24).*15;

disp(GMST)