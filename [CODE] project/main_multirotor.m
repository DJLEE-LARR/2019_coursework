%% main_multirotor.m

clc;
clear;
close all;

m = 1.2;
g = 9.81;

full_DDP = false;
xf = [1;1;1;0;0;0;zeros(6,1)];
% xf = [0;0;0;0;0;0;zeros(6,1)];
udes = [m*g;0;0;0];
h = 0.05;           % time step (sec)

DYNCST = @(x,u,i) quad_dyn_cst(x,u,xf,udes,h,full_DDP);
Op.lims  = [0   25;         % thrust limits (N)
           -10  10;         % torque 1 (Nm)
           -10  10;         % torque 2 (Nm)
           -10  10];        % torque 3 (Nm) 
% Op.lims = [];
Op.plot = -1;               % plot the derivatives as well
Op.parallel = false;
iterMax = 100;
iter = 1;

Tcount = zeros(1,iterMax);
xMPC = zeros(size(xf,1),iterMax);
uMPC = zeros(size(udes,1),iterMax);

% iLQR method
T = iterMax;
% x0      = [0;0;0;0;0;0;zeros(6,1)];   % initial state
x0 = [0;0;0;0;0;0;0;pi;0;zeros(3,1)];
u0      = [m*g*ones(1,T);zeros(3,T)]; % initial controls
[xN_iLQR,u_iLQR]= iLQG_v1(DYNCST, x0, u0, Op);

x_iLQR(:,1) = x0;
for i=1:size(u_iLQR,2)
    x_iLQR(:,i+1) = nominal_dynamics(x_iLQR(:,i),u_iLQR(:,i),h); 
%     x_iLQR(:,i+1) = real_dynamics(x_iLQR(:,i),u_iLQR(:,i),h); 
end

fprintf('finished');
x = reshape(x_iLQR,size(x_iLQR,1),[]);
u = reshape(u_iLQR,size(u_iLQR,1),[]);

Hzn = T;


figure()
sgtitle('iLQR');
subplot(2,2,1)
title('state - x,y,z');
plot(h*(1:Hzn),x(1:3,1:end-1));
xlim([0 h*Hzn]);
legend('x','y','z');
xlabel('time (sec)');
ylabel('Position (m)');

subplot(2,2,2)
title('state - \phi,\theta,\psi');
plot(h*(1:Hzn),RadToDeg(x(7:9,1:end-1)));
xlim([0 h*Hzn]);
legend('\phi','\theta','\psi');
xlabel('time (sec)');
ylabel('Attitude (deg)');

subplot(2,2,3)
title('input');
plot(h*(1:Hzn),u(1:4,:));
xlim([0 h*Hzn]);
legend('f','\tau_x','\tau_y','\tau_z');
xlabel('time (sec)');
ylabel('Input');


% figure()
% sgtitle('iLQR');
% subplot(2,2,1)
% title('state - x,y,z');
% plot(h*(1:Hzn),x(1:3,1:end));
% hold on;
% plot(h*(1:Hzn),x_iLQR(1:3,1:end-1),'--');
% xlim([0 h*Hzn]);
% legend('MPC-x','MPC-y','MPC-z','iLQR-x','iLQR-y','iLQR-z');
% xlabel('time (sec)');
% ylabel('Position (m)');
% 
% subplot(2,2,2)
% title('state - \phi,\theta,\psi');
% plot(h*(1:Hzn),RadToDeg(x(7:9,1:end)));
% hold on;
% plot(h*(1:Hzn),RadToDeg(x_iLQR(7:9,1:end-1)),'--');
% xlim([0 h*Hzn]);
% legend('MPC-\phi','MPC-\theta','MPC-\psi','iLQR-\phi','iLQR-\theta','iLQR-\psi');
% xlabel('time (sec)');
% ylabel('Attitude (deg)');
% 
% subplot(2,2,3)
% title('input');
% plot(h*(1:Hzn),u(1:4,:));
% hold on;
% plot(h*(1:Hzn),u_iLQR(1:4,:),'--');
% xlim([0 h*Hzn]);
% legend('MPC-f','MPC-\tau_x','MPC-\tau_y','MPC-\tau_z','iLQR-f','iLQR-\tau_x','iLQR-\tau_y','iLQR-\tau_z');
% 
% xlabel('time (sec)');
% ylabel('Input');

%% functions
function y = nominal_dynamics(x, u, h)

% === states and controls:
% x = [xq yq zq xqdot yqdot zqdot phi theta psi phidot thetadot psidot].'
% u = [f;tx;ty;tz].' = [thrust;torque].': 4X1

% constants
m = 1.2;       % mass
Jx = 0.02;     % MOI x-axis
Jy = 0.02;     % MOI y-axis 
Jz = 0.02;     % MOI z-axis
g_vec = [0;0;-9.81];
% h  = 0.01;     % h = timestep (seconds)

% controls
f   = u(1,:); % f = thrust
tx  = u(2,:); % tx = x-axis torque : 
ty  = u(3,:); % ty = y-axis torque
tz  = u(4,:); % tz = z-axis torque

F = [zeros(1,size(f,2));zeros(1,size(f,2));f/m];
T = [tx/Jx;ty/Jy;tz/Jz];

phi = x(7,:);     % phi     = roll
theta = x(8,:);   % theta = pitch
psi = x(9,:);     % psi = yaw

xqdot = x(4,:);
yqdot = x(5,:);
zqdot = x(6,:);

phidot = x(10,:);
thetadot = x(11,:);
psidot = x(12,:);

dy = zeros(size(x,1),size(f,2));
for i=1:size(f,2)
    dy(:,i) = [h*[xqdot(:,i);yqdot(:,i);zqdot(:,i)]; h*(g_vec+Rot([phi(:,i),theta(:,i),psi(:,i)])*F(:,i));...
        h*[phidot(:,i);thetadot(:,i);psidot(:,i)]; h*T(:,i)];   % change in state
end
y  = x + dy;                % new state

end

function c = quad_cost(x, u, xf,udes)
% cost function for car-parking problem
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: quadratic cost on final states
% lx: quadratic cost on states

% if isnan(u(1,:)), make u = 0;
final = isnan(u(1,:));  % b/c we designed u(1,end) = NaN
u(:,final)  = 0;        % making nan-values in u = 0

cu  = 1e-1*[0.5 1 1 1];       % control cost coefficients

                          % final state cost coefficients
cf  = 1e0*[2 2 5 ...     % [  x          y      z
           1 1 1 ...     %  phi      theta    psi
           .1 .1 .1 ...     % xdot       ydot   zdot
           .1 .1 .1];       % phidot thetadot psidot]
       
% xf  = zeros(12,1);          % desired final state

                          % running cost coefficients
cx  = 1e0*[2 2 5 ...     % [  x          y      z
           1 1 1 ...     %  phi      theta    psi
           .1 .1 .1 ...     % xdot       ydot   zdot
           .1 .1 .1];       % phidot thetadot psidot]
      
% co = ones(1,12);       

xdes  = xf;                 % desired trajectory
% xo = [0 0 0 zeros(1,9)].';

% control cost
lu    = cu*(u-udes).^2;

% final cost
if any(final)
   llf      = cf*(x(:,final)-xf).^2;
   lf       = double(final);
   lf(final)= llf;
else
   lf    = 0;
end

% running cost
lx = cx*(x-xdes).^2;

% obstacle avoidance cost
% lo = exp(-co*(x-xo).^2);
lo = 0;

% total cost
c     = lu + lx + lf + lo;
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = quad_dyn_cst(x,u,xf,udes,h,full_DDP)
% combine car dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = nominal_dynamics(x,u,h);
    c = quad_cost(x,u,xf,udes);
else
    % state and control indices
    ix = 1:12;
    iu = 13:16;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) nominal_dynamics(xu(ix,:),xu(iu,:),h);
    J       = finite_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    % dynamics second derivatives
    if full_DDP
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = reshape(JJ, [4 6 size(J)]);
        JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);    
    else
        [fxx,fxu,fuu] = deal([]);
    end    
    
    % cost first derivatives
    xu_cost = @(xu) quad_cost(xu(ix,:),xu(iu,:),xf,udes);
    J       = squeeze(finite_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    % cost second derivatives
    xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [x; u]);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    
    [f,c] = deal([]);
end
end

function J = finite_difference(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);
end

function output = Rot(q)
   
    % Euler ZYX angle
    phi = q(1);
    theta = q(2);
    psi = q(3);
    
%     syms phi theta psi
    Rx = [1 0 0;
        0 cos(phi) -sin(phi);
        0 sin(phi) cos(phi)];
    Ry = [cos(theta) 0 sin(theta);
        0 1 0;
        -sin(theta) 0 cos(theta)];
    Rz = [cos(psi) -sin(psi) 0;
        sin(psi) cos(psi) 0;
        0 0 1];

    output = Rz*Ry*Rx;
end
    
% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);
end

function Deg = RadToDeg(Rad)
    Deg = 180/pi*Rad;
end