clear all
%close all

%% parameters
for N = 10:10:100 % number of dimension
J = 2; % J < N
A = diag(1:N);
%A = rand(N);
A = -A; % coefficents for ODE
c = ones(N,1); % initial conditions
g = ones(1,N); %[1,0.5,2,1,1,3,0,0.25]; % weight vector for R^N to compute QoI

%% time
TT = 10;   % total time
dt = 0.1; % time step
tt = 0:dt:TT;  % time vector
n = length(tt)-1; % number of time steps
m = 10; % m < n number of past values to keep
E = expm(A*dt); % evolution matrix

tic

%% projection
P = zeros(N,J);
%P = rand(N,J);
Q = zeros(N,N);
%P = P0;
for i = 1:J
%    P(:,i) = rand(N,1);
    P(i,i) = 1;
end
P(:,J)=ones(N,1);

P1 = P; % Gram-Schmidt(all the vectors in P1 is orthogonal to each other)
P1(:,1) = P(:,1)/norm(P(:,1));
for i = 2:J
    for j = 1:i-1
        p = P1(:,j)/norm(P1(:,j));
        P1(:,i) = P1(:,i) - (P(:,i)'*p)*p;
    end
    P1(:,i) = P1(:,i)/norm(P1(:,i));
end

P = P1(:,1)*P1(:,1)'; % projection matrix
for i = 2:J
    P = P + P1(:,i)*P1(:,i)';
end
P = P';
Q = eye(N) - P; % complementary matrix, orthogonal complement

%% redefine g as the coefficients in the projected space
%g=P1';

%% precomputed evolution matrices
EQ = Q*E;
EP = P*E;
E1 = P1'*E*P1;
% for k = 1:n
%     F(:,:,k) = g*(EQ)^(k);
% end
for k = 1:m
    Fp(:,:,k) = g*(EQ)^(k)*P1;
    E2(:,:,k) = P1'*E*EQ^k*P1;
end
gp = g*P1;

toc
t(N/10) = toc;

tic
%% Initialisation and time Loop
%a = zeros(N,n+1);
cp = zeros(J,n+1);
%ap = zeros(N,n+1);
%Gn = zeros(size(g,1),n+1);
%Mn = zeros(size(g,1),n+1);
Mnp = zeros(size(g,1),n+1);
%noise = zeros(N,n+1);
noise1 = zeros(N,n+1);
%a(:,1) = P*c;
%ap(:,1) = P*c;
cp(:,1) = P1'*c;
b1 = Q*c;
%Gn(:,1) = g*c;
%Mn(:,1) = g*c;
Mnp(:,1) = g*c;
for i = 1:n
% %% exact formula for a
%     a(:,i+1) = EP*a(:,i);
%     for k = 2:i
%         a(:,i+1) = a(:,i+1) + EP*(EQ^(k-1))*a(:,i-k+1);
%     end
%     noise(:,i+1) = EP*EQ^(i-1)*b1;
%     a(:,i+1) = a(:,i+1) + noise(:,i+1);
% %% evolution of G
%     Gn(:,i+1) = g*a(:,i+1); %formula for projection
%     for k = 1:i
%         Gn(:,i+1) = Gn(:,i+1) + F(:,:,k)*a(:,i-k+1);
%     end
%     noise1(:,i+1) = F(:,:,i)*b1;
%     Gn(:,i+1) = Gn(:,i+1) + noise1(:,i+1);
%% approximated formula with limited history
%    ap(:,i+1) = EP*ap(:,i);
    cp(:,i+1) = E1*cp(:,i);
    for k = 2:min(i,m)
%        ap(:,i+1) = ap(:,i+1) + EP*(EQ^(k-1))*ap(:,i-k+1);
        cp(:,i+1) = cp(:,i+1) + E2(:,:,k-1)*cp(:,i-k+1);
    end
%    Mn(:,i+1) = g*ap(:,i+1); %formula for projection
    Mnp(:,i+1) = gp*cp(:,i+1); %formula for projection
    for k = 1:min(i,m)
%        Mn(:,i+1) = Mn(:,i+1) + F(k)*ap(:,i-k+1);
        Mnp(:,i+1) = Mnp(:,i+1) + Fp(:,:,k)*cp(:,i-k+1);
    end
end

toc
t1(N/10) = toc;
%% Direct method
%odefun = @(t,x) A*x;
%[t,y]=ode113(odefun, tt ,c); % solve the equation
%G = g*y(end,:)'; % Gn = g*cn
%G0 = g*y';

%% analytical solution
Gexact = g*(c.*exp(diag(A)*tt));

% gg = E^n*g'; % gn = E^n*g
% G1 = gg'*c; % Gn = gn*c0

loglog(tt,Gexact,'y-o','Displayname', 'Analytical solution');
hold on;
%loglog(tt,Gn,'b-o','Displayname', 'Full model');
%loglog(tt,G0,'k-o','Displayname', 'ODE solution'); 
loglog(tt,Mnp,'r-o','Displayname', 'Reduced model');

err(N/10)=norm(Gexact-Mnp)/norm(Gexact);
end

figure;
i = 10:10:100;
plot(i,t,'-o')
hold on 
plot(i,t1,'-o')
xlabel('N');
ylabel('time');
grid on
title('time used for precomputed and approximate for different N');
legend('precomputed', 'approximate');

figure;
plot(i,err,'-o')
xlabel('N');
ylabel('error');
grid on
title('error between approximation and exact value');

