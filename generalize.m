clear all

%% parameters
N = 8; % number of dimension
J = 2; % J < N
A = diag(1:N); 
A = -A; % coefficents for ODE
c = rand(N,1); % initial conditions
g = [1,0.5,2,1,1,3,0,0.25]; % weight vector for R^N to compute QoI

%% time
TT = 10;
dt = 0.1;
tt = 0:dt:TT;
n = length(tt)-1; % times of evolution
m = 1; % m < n number of past values to keep
E = expm(A*dt); % evolution matrix

%% projection
for i = 1:J
   P(:,i) = rand(N,1);
%    P(:,i) = zeros(N,1);
  %  P(i,i) = 1;
end

P1 = P; % Gram£­Schmidt(all the vectors in P1 is orthogonal to each other)
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


%% exact formula
gnp = g*( P*(E^n) )*c;
for k = 1:n
    gnp = gnp + g*(((Q*E)^(k-1))*Q*E*P*(E^(n-k)) )*c; %formula for projection
end
gnq = g*(Q*E)^n*Q*c; %orthogonal part, need the information of c from initial time
gn = gnp + gnq;
Gn = gn %*c

%% approximated formula with limited history
g1n = g*P*(E^n)*c;
for k = 1:m
    g1n = g1n + g*(((Q*E)^(k-1))*Q*E*P*(E^(n-k)) )*c; %formula for projection, E^(n-k)*c = c_(n-k)(markov) we only need to have the information of c at n-k
end
Mn = g1n %the approximated value


%% Direct method
odefun = @(t,x) A*x;
[t,y]=ode113(odefun, tt ,c); % solve the equation
G = g*y(end,:)' % Gn = g*cn

gg = E^n*g'; % gn = E^n*g
G1 = gg'*c % Gn = gn*c0

%% analytical solution
Gexact = g*(c.*exp(diag(A)*TT))

disp('Errors:')
disp([1-Gn/Gexact 1-Mn/Gexact 1-G/Gexact 1-G1/Gexact] )

