% - Solves FitzHugh Nagumo travelling wave equation
% - Requires optimization toolbox and external routines FHNeqn_w_periodic.m
% and fourdif.m

% Copyright (C) 2023 Paul Carter
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <https://www.gnu.org/licenses/>.

%% setup

par.N = 1000;      % number of mesh points, must be even! 
par.L = 1000;       % domain truncation

%  spatial grid
N = par.N;
L = par.L;
par.h = L/(N-1); h = par.h;
x = (0:N-1)'*h;

%% initial conditions & parameters

% params1 [§8.1 Figure 6]
 par.a = 0.0997212;
 par.eps =  0.0021398784312;
 par.gamma = 3.5;
 load('initial_pulse_1');
 c = 0.548;
 lambdac0 = -0.0194;
 par.eta = 0.0;

% % params2 [§8.2.1 Figure 7a]
% par.a = 0.167095;
% par.eps =  0.0021398784312;
% par.gamma = 0.5;
% load('initial_pulse_2');
% c = 0.4446;
% lambdac0 = -0.0408;
% par.eta = 0.1; % exponential weight for eigenvalue problem

% % params3 [§8.2.2 Figure 7b]
% par.a = 0.0058792294665;
% par.eps =  0.0021398784312;
% par.gamma = 0.5;
% load('initial_pulse_3');
% c = 0.6864;
% lambdac0 = -0.0374;
% par.eta = 0.1; % exponential weight for eigenvalue problem

U  = Uout; 
U0 = U;

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',100000,'Algorithm','levenberg-marquardt');

% call solve
[Uout,fval] = fsolve(@(U) FHNeqn_w_periodic(U,U0,par),U,options);

%% plot results

cout = Uout(end)  % wave speed

figure(1);
plot(x,Uout(1:N));
xlabel('x');
ylabel('u');
title('Pulse on the line');
drawnow;
 
%% compute leading eigenfunction

[F,J,Jw] = FHNeqn_w_periodic(Uout,Uout,par);
[V1,DD1] = eigs(Jw,100,1); 
evalsn = diag(DD1);

% manually extract critical eigenvalue/eigenfunction
eval_loc = find(abs(evalsn-lambdac0)<0.001);
lambdac = evalsn(eval_loc) % critical eigenvalue
v1 = V1(:,eval_loc);
u1 = real(v1(1:N));
figure(2);
plot(x,Uout(1:N));hold on;
plot(x,max(Uout(1:N))/max(abs(u1))*u1); hold off;
xlabel('x');
ylabel('u');
title('Critical eigenfunction in weighted space');
