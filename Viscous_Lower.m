%--------------------------------------------------------------------------
% Bounding the minimum wave speed of the viscous FKKP-Burgers equation: 
%
%    T_t - T_xx + (u*T)_x = T*(1-T),
%   u_t - nu*u_xx + u*u_x = rho*T*(1-T),    rho,nu > 0
%
% This code produces a lower bound on the (conjectured) minimum wave speed 
% and is associated to the paper "The speed of traveling waves in a 
% FKPP-Burgers system" by Jason J. Bramburger and Christopher Henderson 
% (2020). This script is used to generate the data in Table 1.
%--------------------------------------------------------------------------

% Access YALMIP and mosek directories
addpath(genpath('YALMIP-master')) 
addpath(genpath('mosek')) 

% Clean workspace
clear all
close all
clc

format long

% Differential equation Parameters
nu = 1;
rho = 10^6; 

% Bounding method parameters
d = 10; %Degree of H(T,U,V)
lambda = 1;
    
%Bisection Method
cleft = 2;
cright = 1000;

while abs(cright - cleft) >= 1e-4

    cmid = 0.5*(cleft + cright);
    flag = lowerbnd(cmid,rho,nu,d,lambda);

    if flag == 0
        cleft = cmid; 
    else
       cright = cmid;
    end
end

%Print Results
fprintf('The lower bound on the minimum speed for rho = %d and nu = %d is %f found using degree %d polynomials.\n',rho,nu,cmid,d)

%%
function flag = lowerbnd(c,rho,nu,d,lambda)

% SDP variables
sdpvar T U V

% Source fixed point value
u0 = c + rho - sqrt(c^2 + rho^2);
v0 = c - u0;

% Epsilon value
eps = 1e-4;

% Auxiliary function
[H, cH] = polynomial([T U V], d);

% S Procedure
d2 = d + 4;
[s1, c1] = polynomial([T U V], d2);
[s2, c2] = polynomial([T U V], d2);
[s3, c3] = polynomial([T U V], d2);
[s4, c4] = polynomial([T U V], d2);
[s5, c5] = polynomial([T U V], d2);
[s6, c6] = polynomial([T U V], d2);

% Derivatives
dHdT = jacobian(H,T);
dHdU = jacobian(H,U);
dHdV = jacobian(H,V);

%Constraints
cons = [];
cons = [cons, replace(H, [T U V], [0 0 0]) == 0]; 
cons = [cons, replace(H, [T U V], [1 1 1]) <= -eps];
cons = [cons, sos(-lambda*(dHdT*(-c*T + U*u0*T + v0*V) + dHdU*((-c*U*u0 + 0.5*(u0^2)*U^2 + rho*v0*V)/(nu*u0)) + dHdV*(T*(T-1)/v0)) - H - T*(1-T)*s1 - U*(1-U)*s2 - V*(1-V)*s3)];
cons = [cons, sos(H - T*(0.05-T)*s4 - U*(0.05-U)*s5 - V*(0.05-V)*s6)];
cons = [cons, sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6)];

%SOS Solver
ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
sol = solvesos(cons,[],ops,[cH; c1; c2; c3; c4; c5; c6])

%Return whether solvesos failed or succeeded
flag = sol.problem;

end











