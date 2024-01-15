%Orbital Mechanics Simulator For Primary Satellite and Orbial Debris given
%2LE Data

%Patrick Williams
%Aerospace and Mechanical Engineering MS Student
%Oklahoma State University
%patrick.j.williams@okstate.edu
%Autonomous Physics Group autophysics.net

%Last Updated: Jan 15, 2024


%Solving for the initial conditions from 
%Analytical Mechanics of Space Systems 4th Ed.
%Chapter 9
function [T, X0] = Sat_X0(OE,c)
        
        Eps = 2.22044604925031e-8;                 % Machine epsilon

        %OE = [a,e,i,RAAN,omega,M]

        a = OE(1)*1000; %km to m
        e = OE(2);
        i = deg2rad(OE(3));
        RAAN = deg2rad(OE(4));
        Omega = RAAN;
        omega = deg2rad(OE(5));
        w = omega;
        M = deg2rad(OE(6));

        mu = c.G * (c.M_E);
        



        E   = SolveKepler(M,e,Eps);                % Eccentric Anomaly from Mean Anomaly
        M   = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));  % True Anomaly obtained from Eccentric Anomaly
        p = a*(1-e^2); %m


        %Orbit Period
        T       = 2*pi*sqrt(a^3/mu); % s
        %tfin    = 1*T;
        f = M;%E0;
        %r = p/(1+e*cosd(f));
        %theta = w + f;


        x0   = (p*cos(M)) / (1 + e*cos(M));         % x-coordinate (m)
        y0   = (p*sin(M)) / (1 + e*cos(M));         % y-coordinate (m)
        z0   = 0;                                   % z-coordinate (m)
        u0  = -(mu/p)^(1/2) * sin(M);              % velocity in x (m/s)
        v0  = (mu/p)^(1/2) * (e + cos(M));         % velocity in y (m/s)
        w0  = 0; 

        %3-1-3 Euler Angle Rotation Matrix
        t_rot = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i) ...
            (-1)*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i) ...
            sin(Omega)*sin(i); ...
            sin(Omega)*cos(w)+cos(Omega)*sin(w )*cos(i ) ...
            (-1)*sin(Omega )*sin(w )+cos(Omega )*cos(w )*cos(i ) ...
            (-1)*cos(Omega )*sin(i ); ...
            sin(w )*sin(i )  cos(w )*sin(i )  cos(i )];

        R = t_rot*[x0 y0 z0]';
        V = t_rot*[u0 v0 w0]';

        X0 = [R(1), R(2), R(3), V(1), V(2), V(3)];
end