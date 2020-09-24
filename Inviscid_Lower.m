%--------------------------------------------------------------------------
% Bounding the minimum wave speed of the inviscid FKKP-Burgers equation: 
%
%    T_t - T_xx + (u*T)_x = T*(1-T),
%             u_t + u*u_x = rho*T*(1-T),    rho > 0
%
% This code produces a lower bound on the minimum wave speed and is 
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
cleft = 110;
cright = 120;

while abs(cright - cleft) >= 1e-4

    cmid = 0.5*(cleft + cright);
    flag = lowerbnd(cmid,rho,d,lambda);

    if flag == 0
        cleft = cmid; 
    else
       cright = cmid;
    end
end

%Print Results
fprintf('The lower bound on the minimum speed for rho = %d is %f found using degree %d polynomials.\n',rho,cmid,d)

%%
function flag = lowerbnd(c,rho,d,lambda)

% SDP variables
sdpvar T U

% Epsilon value
eps = 1e-4;

% Source fixed point value
u0 = c + rho - sqrt(c^2 + rho^2);

% Auxiliary function
[H, cH] = polynomial([T U],d);

% S Procedure
d2 = d + 4;
[s1, c1] = polynomial([T U], d2);
[s2, c2] = polynomial([T U], d2);
[s4, c4] = polynomial(U, d2);

% Derivatives
dHdT = jacobian(H,T);
dHdU = jacobian(H,U);

% Function replacements
HT1 = replace(H,T,1); %H(0,U)

%Constraints
cons = [];
cons = [cons, replace(H,[T U], [0 0]) >= eps]; 
cons = [cons, replace(H,[T U], [1 1]) == 0];
cons = [cons, sos(-lambda*(dHdT*(-c*T + U*u0*T + U*u0*(2*c-U*u0)/(2*rho))*(c-U*u0) + dHdU*(-rho*T*(1-T))/u0) + H*(c-U*u0) - T*(1-T)*s1 - U*(1-U)*s2)];
cons = [cons, sos(-HT1 - eps*U*(1-U) - U*(1-U)*s4)];
cons = [cons, sos(s1), sos(s2), sos(s4)];

%SOS Solver
ops = sdpsettings('solver','mosek','verbose',0);
sol = solvesos(cons,[],ops,[cH; c1; c2; c4])

%Return whether solvesos failed or succeeded
flag = sol.problem;

end











