%clear all

%% parameters1
N = 10; % number of dimension
J = 2; % J < N
A = diag(1:N); 
A = -A; % coefficents for ODE
odefun = @(t,x) A*x;
c0 = ones(N,1); % initial conditions
g = ones(1,N); % weight vector for R^N to compute QoI

%% time
TT = 1;
dt = 0.1;
tt = 0:dt:TT;
n = length(tt)-1; % times of evolution
m = 5; % m < n number of past values to keep

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

%% precomputed evolution matrices

E1(:,:) = P1' * evolreduced(odefun,P1,dt); %once evolution


for i = 1:m
    E2(:,:,i) = P1' * (evolreduced(odefun,Q,dt))^i * evolreduced(odefun,P1,dt);
    Fp(:,:,i) = g * (Q * (evolreduced(odefun,Q,dt))^(i-1) * evolreduced(odefun,P1,dt));
end

  
%% Initialisation and time Loop
cp = zeros(J,n+1);
cp(:,1) = P1'*c0;
Mnp = zeros(size(g,1),n+1);
Mnp(:,1) = g*c0;

%% Approximate value
for i = 1:n
        cp(:,i+1) = E1*cp(:,i);
    for k = 2:min(i,m)
        cp(:,i+1) = cp(:,i+1) + E2(:,:,k-1)*cp(:,i-k+1);
    end
    Mnp(:,i+1) = g*P1*cp(:,i+1); %formula for projection
    for k = 1:min(i,m)
        Mnp(:,i+1) = Mnp(:,i+1) + Fp(:,:,k)*cp(:,i-k+1);
    end
end


