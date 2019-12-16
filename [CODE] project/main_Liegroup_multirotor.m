%% main_multirotor.m

clc;
clear;
close all;
%% Assumption
% 1) All states used in function inputs and outputs are vectors.
%     -> SE(3) components are all mapped to a vector with inverse-cayley map.
% 2) States should be mapped back into original SE(3) with cayley map for
% computation.
% 3) States should be remapped back into a vector with inverse-cayley map.

%%
m = 1.2;
g = 9.81;

% state
    % x = [g,V] \in SE(3) X R^6
    % xc = cay(x) \in R^12
    % g = [R p; 0 1], V = [w v]
full_DDP = false;
xf.g = [eye(3) [1;1;1];
        zeros(1,3) 1];
xf.V = zeros(6,1);
xfcay = [cay_inv(xf.g);xf.V];

udes = [0;0;0;m*g];
h = 0.05;           % time step (sec)

DYNCST = @(xcay,u,i) quad_dyn_cst(xcay,u,xfcay,udes,h,full_DDP);
Op.lims  = [-10  10;         % torque 1 (Nm)
            -10  10;         % torque 2 (Nm)
            -10  10;         % torque 3 (Nm) 
             0   25];        % thrust limits (N)
       
% Op.lims = [];
% Op.plot = -1;               % plot the derivatives as well

iterMax = 100;
iter = 1;

Tcount = zeros(1,iterMax);


% iLQR method
T = iterMax;
x0.g = [eul2rotm([0 0 0]) [0;0;0];
        zeros(1,3) 1];
x0.V = zeros(6,1);
x0cay = [cay_inv(x0.g);x0.V];
u0    = [zeros(3,T);m*g*ones(1,T)]; % initial controls
[xN_iLQR,u_iLQR]= iLQG_v1(DYNCST, x0cay, u0, Op);

% x_iLQR(:,1) = x0cay;
% for i=1:size(u_iLQR,2)
%     x_iLQR(:,i+1) = nominal_dynamics(x_iLQR(:,i),u_iLQR(:,i),h); 
% %     x_iLQR(:,i+1) = real_dynamics(x_iLQR(:,i),u_iLQR(:,i),h); 
% end

fprintf('finished');
Hzn = T;

% cayley map
for i=1:size(xN_iLQR,2)
    SE = cay(xN_iLQR(1:6,i));
    x_cay(1:3,i) = SE(1:3,4);               % position
    x_cay(4:6,i) = rotm2eul(SE(1:3,1:3),'ZYX'); % euler angle
    x_cay(7:12,i) = xN_iLQR(7:12,i);         % generalized velocity
end

figure()
sgtitle('iLQR in Lie Group');
subplot(2,2,1)
title('state - x,y,z');
plot(h*(1:Hzn),x_cay(1:3,1:end-1),'-');
xlim([0 h*Hzn]);
legend('x','y','z');
xlabel('time (sec)');
ylabel('Position (m)');

subplot(2,2,2)
title('state - \phi,\theta,\psi');
plot(h*(1:Hzn),RadToDeg(x_cay(4:6,1:end-1)),'-');
xlim([0 h*Hzn]);
legend('\phi','\theta','\psi');
xlabel('time (sec)');
ylabel('Attitude (deg)');

subplot(2,2,3)
title('input');
plot(h*(1:Hzn),u_iLQR(1:4,:),'-');
xlim([0 h*Hzn]);
legend('\tau_x','\tau_y','\tau_z','f');
xlabel('time (sec)');
ylabel('Input');

%% functions
function out = nominal_dynamics(xcay, u, h)

% === states and controls:
% x = [g V].' \in SE(3) X R^6
    % g = [R p; 0 1], V = [v w]
    % xcay = [cay_inv(x.g);x.V] \in R^12
% u = [tx;ty;tz;f].' = [torque;thrust].': 4X1

% constants
m = 1.2;       % mass
Jx = 0.02;     % MOI x-axis
Jy = 0.02;     % MOI y-axis 
Jz = 0.02;     % MOI z-axis
g_vec = [0;0;9.81];
J = diag([Jx,Jy,Jz]);
M = diag([Jx,Jy,Jz,m,m,m]);
B = [eye(3) zeros(3,1);
     zeros(3,3) [0;0;1]];
    
% what if there exists more than 1 column of x?
sz = size(xcay,2);
out_temp = zeros(12,sz);

for idx=1:sz
    g = cay(xcay(1:6,idx));
    R = g(1:3,1:3);
    % p = g(1:3,4);
    w = xcay(7:9,idx);
    v = xcay(10:12,idx);
    
    b = [cross(w,J*w);cross(m*w,v) + m*R.'*g_vec];
   
    a_k = inv(M)*(-b + B*u(:,idx));
    f_k = h*[[w;v] + h*a_k;a_k];

    % output (next step state)
    xnew.g = g*cay(f_k(1:6,1));
    xnew.V = [w;v] + f_k(7:12,1);

    out_temp(:,idx) = [cay_inv(xnew.g);xnew.V];
end
    out = out_temp;
end

function c = quad_cost(xcay, u, xfcay, udes)
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: quadratic cost on final states
% lx: quadratic cost on states

% if isnan(u(1,:)), make u = 0;
final = isnan(u(1,:));  % b/c we designed u(1,end) = NaN
u(:,final)  = 0;        % making nan-values in u = 0

% Ru: control cost coefficients
cu  = 1e0*[1 1 1 0.5];       

% Qf: final state cost coefficients                         
cf  = 1e1*[5 5 5 ...     % [ cay-1(q1^-1q2)
           5 5 5 ...     %   V];
           .1 .1 .1 ...     
           .1 .1 .1];       
       
% xf  = zeros(12,1);          % desired final state

% Qv: running cost coefficients
cx  = 1e1*[5 5 5 ...     % [ cay-1(q1^-1q2)
           5 5 5 ...     %   V];
           .1 .1 .1 ...     
           .1 .1 .1];       
      
% co = ones(1,12);       

xdes  = xfcay;                 % desired trajectory
% xo = [0 0 0 zeros(1,9)].';

% control cost
lu    = cu*(u-udes).^2;

% final cost
if any(final)
   xf.g = cay(xfcay(1:6,1));
   xf.V = xfcay(7:12,1);
   xend.g = cay(xcay(1:6,final));
   xend.V = xcay(7:12,final);
   
   err = [cay_inv(inv(xf.g)*xend.g);xend.V-xf.V];
   llf      = 1/2*cf*(err).^2;
   lf       = double(final);
   lf(final)= llf;
else
   lf    = 0;
end

% running cost
xs = size(xcay,2);
for idx = 1:xs
xd_.g = cay(xdes(1:6,1));
xd_.V = xdes(7:12,1);
x_.g = cay(xcay(1:6,idx));
x_.V = xcay(7:12,idx);
   
err_ = [cay_inv(inv(xd_.g)*x_.g);x_.V-xd_.V];

lx(:,idx) = cx*(err_).^2;
end

% total cost
c     = lu + lx + lf;
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = quad_dyn_cst(xcay,u,xfcay,udes,h,full_DDP)
% combine car dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = nominal_dynamics(xcay,u,h);
    c = quad_cost(xcay,u,xfcay,udes);
else
    % state and control indices
    ix = 1:12;
    iu = 13:16;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) nominal_dynamics(xu(ix,:),xu(iu,:),h);
    J       = finite_difference_dynamics(xu_dyn, [xcay; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    % dynamics second derivatives
    if full_DDP
        xu_Jcst = @(xu) finite_difference_dynamics(xu_dyn, xu);
        JJ      = finite_difference_dynamics(xu_Jcst, [xcay; u]);
        JJ      = reshape(JJ, [4 6 size(J)]);
        JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);    
    else
        [fxx,fxu,fuu] = deal([]);
    end    
    
    Qq = diag([1 1 1 1 1 1]);
    Qv = diag([1 1 1 1 1 1]);
    % cost first derivatives
    xu_cost = @(xu) quad_cost(xu(ix,:),xu(iu,:),xfcay,udes);
    J       = squeeze(finite_difference(xu_cost, [xcay; u]));
%     cx      = J(ix,:);
    cu      = J(iu,:);
    sx = size(xcay,2);    
    for i=1:sx
        Delta_q = Deltaq(cay(xfcay(1:6,1)),cay(xcay(1:6,i)));
        cx(:,i) = [dcay_inv(-Delta_q).'*Qq*Delta_q;
                  Qv*(xcay(7:12,i)-xfcay(7:12,1))];
    end
    % cost second derivatives
    xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [xcay; u]);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
%     cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    for i=1:sx
        Delta_q = Deltaq(cay(xfcay(1:6,1)),cay(xcay(1:6,i)));
        cxx(:,:,i) = [dcay_inv(-Delta_q).'*Qq*dcay_inv(-Delta_q) zeros(6,6);
                  zeros(6,6) Qv];
    end
    [f,c] = deal([]);
end
end

function J = finite_difference_dynamics(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end
% n = 16, K = 101
[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
% only changing the line below would be fine
    % minus between cayley-inverse mapped states cannot be done directly
    for idx1 = 1:K
        g1 = cay(Y(1:6,idx1,1));
        R1 = g1(1:3,1:3); p1 = g1(1:3,4);
        inv_g1 = [R1.' -R1.'*p1;
                  zeros(1,3) 1];
        for idx2 = 2:n+1
            g2 = cay(Y(1:6,idx1,idx2));
            J(1:6,idx1,idx2-1) = cay_inv(inv_g1*g2);
        end
    end    
J(7:12,:,:) = pp(Y(7:end,:,2:end), -Y(7:end,:,1)) / h;

J       = permute(J, [1 3 2]);
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

function out = cay(V)
% if V \in R^6
w = V(1:3,1);
v = V(4:6,1);
out = [eye(3)+4/(4+norm(w,2)^2)*(hat(w) + hat(w)^2/2) 4/(4+norm(w,2)^2)*(2*eye(3) + hat(w))*v;
       zeros(1,3) 1];
end

function out = dcay(V)
w = V(1:3,1);
v = V(4:6,1);

out(1:3,1:3) = 2/(4+norm(w,2)^2)*(2*eye(3)+hat(w));
out(1:3,4:6) = zeros(3,3);
out(4:6,1:3) = 1/(4+norm(w,2)^2)*hat(v)*(2*eye(3)+hat(w));
out(4:6,4:6) = eye(3) + 1/(4+norm(w,2)^2)*(2*hat(w) + hat(w)^2);
end

function out = dcay_inv(V)
w = V(1:3,1);
v = V(4:6,1);

out(1:3,1:3) = eye(3) - 1/2*hat(w) + 1/4*w*w.';
out(1:3,4:6) = zeros(3,3);
out(4:6,1:3) = -1/2*(eye(3) - 1/2*hat(w))*hat(v);
out(4:6,4:6) = eye(3) - 1/2*hat(w);
end

function out = Deltaq(q1,q2)
% cay_inv(q1^-1q2)
R1 = q1(1:3,1:3);
p1 = q1(1:3,4);
inv_q1 = [R1.' -R1.'*p1;
          zeros(1,3) 1];
out = cay_inv(inv_q1*q2);
end

function out = cay_inv(g)
    R = g(1:3,1:3);
    p = g(1:3,4);
out = [vee(-2*inv((eye(3)+R))*(eye(3)-R));
       inv((eye(3)+R))*p];
end

function out = hat(R3)
   out = [0 -R3(3) R3(2);
          R3(3) 0 -R3(1);
         -R3(2) R3(1) 0];
end

function out = vee(skew)
    % skew = [0 -a3 a2;
%            a3 0 -a1;
%           -a2 a1 0]
   out = [skew(3,2);
          skew(1,3);
          skew(2,1)];
end