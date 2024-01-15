function E = SolveKepler(M,e,Eps)
% Iterative Solution of Kepler's equation (M = E-e*sin(E))
%==========================================================================
% Inputs:    M   =  Mean Anomaly (rad)
%            e   =  Eccentricity
%            Eps =  Machine epsilon
% Output:    E   =  Eccentric Anomaly (rad)
%==========================================================================
E0  = M;
E1 = E0 - (E0 - e*sin(E0) - M)/(1  -e*cos(E0));
while (abs(E1-E0) > Eps)
    E0 = E1;
    E1 = E0 - (E0 - e*sin(E0) - M)/(1 - e*cos(E0));
end
E = E1;
end