function [F,J,Jw] = FHNeqn_w_periodic(U,U0,par)
% returns the right-hand side and its Jacobian of FitzHugh-Nagumo equation
% 0 = - cu_x + u_xx + u(1-u)(u-a) - w
% 0 = - cw_x + epsilon(u-gamma w)

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

%% parameters

a = par.a;
eps = par.eps;
eta = par.eta;
gamma = par.gamma;
N = par.N;
h = par.h;
L = par.L;

u = U(1:N);
w = U(N+1:2*N);
u0 = U0(1:N);
w0 = U0(N+1:2*N);
c = U(end);

%% compute linear FitzHugh-Nagumo operator & differentiation matrices

z = sparse(N,N);
e = speye(N,N); % identity matrix

% spectral differentiation matrices
D = fourdif(N,1)*2*pi/L; % d_x
D2 = fourdif(N,2)*(2*pi/L)^2; % d_xx

LD = [D z; z D];
LN = -c*LD + [D2 z; z z] + [z -e; eps*e -eps*gamma*e];

%% FHN right hand side
F = [LN*[u;w] + [-u.^3+(1+a)*u.^2-a*u; sparse(N,1)];...
    (LD*[u;w])'*([u;w] - [u0;w0])];

%% Jacobian
if nargout > 1
    J = [LN + spdiags([-3*u.^2+2*(1+a)*u-a; sparse(N,1)],0,2*N,2*N), -LD*[u;w];...
        (LD*[u;w])'+ [u;w]'*LD - [u0;w0]'*LD, 0];
    if nargout>2
        % weighted operators
        LDw = LD+eta*speye(2*N);
        D2w = D2+2*eta*D+eta^2*speye(N);
        LNw = -c*LDw + [D2w z; z z] + [z -e; eps*e -eps*gamma*e];
        Jw = [LNw + spdiags([-3*u.^2+2*(1+a)*u-a*ones(N,1); sparse(N,1)],0,2*N,2*N)]; % Jacobian for weighted problem
    end
end

end
