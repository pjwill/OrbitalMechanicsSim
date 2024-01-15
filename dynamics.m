%Dynamics Calculations For Orbital Debris Simulator

%Patrick Williams
%Aerospace and Mechanical Engineering MS Student
%Oklahoma State University
%patrick.j.williams@okstate.edu
%Autonomous Physics Group autophysics.net

%Last Updated: Jan 15, 2024


function [dX] = dynamics(~,X0,c)
rx = X0(1);
ry = X0(2);
rz = X0(3);
vx = X0(4);
vy = X0(5);
vz = X0(6);


mu = c.mu;
r  = sqrt(rx.^2+ry.^2+rz.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Update velocities based on orbital eom%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX(1,1) = vx; %m/s
dX(2,1) = vy; %m/s
dX(3,1) = vz; %m/s

dX(4,1) = (-mu/(r^3).*rx + c.fd_2); %Eqn 9.79
dX(5,1) = (-mu/(r^3).*ry + c.fd_2); %Eqn 9.80
dX(6,1) = (-mu/(r^3).*rz + c.fd_2); %Eqn 9.81



end