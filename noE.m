clear all

%% parameters1
N = 10; % number of dimension
J = 2; % J < N
A = diag(1:N); 
A = -A; % coefficents for ODE
odefun = @(t,x) A*x;
c0 = ones(N,1); % initial conditions
g = ones(N,1)'; % weight vector for R^N to compute QoI

%% time
TT = 1;
dt = 0.1;
tt = 0:dt:TT;
n = length(tt)-1; % times of evolution
m = 5; % m < n number of past values to keep

%% projection
for i = 1:J
%    P(:,i) = rand(N,1);
    P(:,i) = zeros(N,1);
    P(i,i) = 1;
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

%% parameter2

for i = 1:n
c(i,:) = evol(odefun,c0,tt(i+1));
end
c = c';
c = [c0,c];
a = P * c;
b1 = Q * c0;
for i = 1:n
    for j = 1:n
    EP(j,:,i) = evol(odefun,a(:,i),tt(j+1));
    end
end

for i = 1:n
    EQ1(i,:) = evol(odefun,b1,tt(i+1));
end
%% Exact formula
Gn(:,1) = g*c0;
for i = 1:n
    Gn(:,i+1) = g*a(:,i+1);
    for k = 1:i
        Gn(:,i) = Gn(:,i) + g*(Q)^(k)*EP(k,:,i-k+1)';
    end
    noise(:,i+1) = g*(Q)^(k)*EQ1(i,:)';
    Gn(:,i+1) = Gn(:,i+1) + noise(:,i+1);
end


