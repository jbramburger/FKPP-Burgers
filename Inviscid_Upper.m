%--------------------------------------------------------------------------
% Bounding the minimum wave speed of the inviscid FKKP-Burgers equation: 
%
%    T_t - T_xx + (u*T)_x = T*(1-T),
%             u_t + u*u_x = rho*T*(1-T),    rho > 0
%
% This code produces an upper bound on the minimum wave speed and is 
% associated to the paper "The speed of traveling waves in a FKPP-Burgers 
% system" by Jason J. Bramburger and Christopher Henderson (2020). 
% This script is used to generate the data in Table 1.
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
rho = 10^6; 

% Bounding method parameters
d = 10; %Degree of H(T,U)
lambda = 1;

%Bisection Method
cleft = 100;
cright = 120;

while abs(cright - cleft) >= 1e-4

    cmid = 0.5*(cleft + cright);
    flag = upperbnd(cmid,rho,d,lambda);

    if flag == 0
       cright = cmid; 
    else
       cleft = cmid;
    end
end

%Print Results
fprintf('The upper bound on the minimum speed for rho = %d is %f found using degree %d polynomials.\n',rho,cmid,d)

%%
function flag = upperbnd(c,rho,d,lambda)

% SDP variables
sdpvar T U

% Source fixed point value
u0 = c + rho - sqrt(c^2 + rho^2);

% Epsilon value
eps = 1e-4;

% Auxiliary function
[H, cH] = polynomial([T U],d);

% S Procedure
d2 = d + 4;
[s1, c1] = polynomial([T U], d2);
[s2, c2] = polynomial([T U], d2);
[s3, c3] = polynomial(T, d2);

% Derivatives
dHdT = jacobian(H,T);
dHdU = jacobian(H,U);

% Function replacements
HU0 = replace(H,U,0); %H(T,0)

%Constraints
cons = [];
cons = [cons, replace(H,[T U], [0 0]) == 0]; 
cons = [cons, replace(H,[T U], [1 1]) <= -eps];
cons = [cons, sos(-lambda*(dHdT*(-c*T + U*u0*T + U*u0*(2*c-U*u0)/(2*rho))*(c-U*u0) + dHdU*(-rho*T*(1-T))/u0) - H*(c-U*u0) - T*(1-T)*s1 - U*(1-U)*s2)];
cons = [cons, sos(HU0 - eps*T*(1-T) - T*(1-T)*s3)];
cons = [cons, sos(s1), sos(s2), sos(s3)];

%SOS Solver
ops = sdpsettings('solver','mosek','verbose',0);
sol = solvesos(cons,[],ops,[cH; c1; c2; c3])

%Return whether solvesos failed or succeeded
flag = sol.problem;

end











