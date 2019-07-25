%clear all

%% parameters1
N = 100; % number of dimension
J = 1; % J < N
A = diag(1:N); 
A = -A; % coefficents for ODE
odefun = @(t,x) A*(sin(x).*x);
c0 = ones(N,1); % initial conditions
g = ones(N,1)'; % weight vector for R^N to compute QoI

%% time
TT = 10;
dt = 0.1;
tt = 0:dt:TT;
n = length(tt)-1; % times of evolution
m = 1; % m < n number of past values to keep

%% projection
P = zeros(N,J);
%P = rand(N,J);
%Q = zeros(N,N);
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

%P = P1(:,1)*P1(:,1)'; % projection matrix
%for i = 2:J
%    P = P + P1(:,i)*P1(:,i)';
%end
%P = P';
%Q = eye(N) - P; % complementary matrix, orthogonal complement

%% Initialisation and time Loop
cp = zeros(J,n+1);
Gn = zeros(size(g,1),n+1);
cp(:,1) = P1'*c0;
%b1 = Q*c0;
Gn(:,1) = g*c0;
gp = g*P1;

tic
%% approximated formula with limited history and no noise
for i = 1:n
    cp(:,i+1) = evolprojected(odefun,P1,cp(:,i),dt);
    for k = 2:min(i,m)
        cp(:,i+1) = cp(:,i+1) + evolorthogonal(odefun,P1,cp(:,i-k+1),dt,k-1);
    end
    Gn(:,i+1) = gp*cp(:,i+1); %formula for projection
    for k = 1:min(i,m)
        Gn(:,i+1) = Gn(:,i+1) + evolF(odefun,g,P1,cp(:,i-k+1),dt,k);
    end
end
toc

tic
%% Direct method
%odefun = @(t,x) A*x;
[t,y]=ode45(odefun, tt ,c0); % solve the equation
%G = g*y(end,:)'; % Gn = g*cn
Gexact = g*y';
toc

norm(Gn-Gexact)/norm(Gexact)


    
