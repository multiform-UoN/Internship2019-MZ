

%% parameters1
N = 10; % number of dimension
J = 2; % J < N
A = diag(1:N); 
A = -A; % coefficents for ODE
c = rand(N,1); % initial conditions
g = [1,0.5,2,1,1,3,0,0.25]; % weight vector for R^N to compute QoI

%% time
TT = 1;
dt = 0.1;
tt = 0:dt:TT;
n = length(tt)-1; % times of evolution
m = 1; % m < n number of past values to keep
E = expm(A*dt); % evolution matrix

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

%% parameters2
EQ = Q*E;
EP = P*E;
F = @(k) g*(EQ)^(k);

%% exact formula
a(:,1) = P*c;
b1 = Q*c;
for i = 1:n
    a(:,i+1) = EP*a(:,i);
    for k = 2:i
        a(:,i+1) = a(:,i+1) + EP*(EQ^(k-1))*a(:,i-k+1);
    end
    noise(:,i) = EP*EQ^(i-1)*b1;
    a(:,i+1) = a(:,i+1) + noise(:,i);
end


for j = 0:n
    Gn1(j+1) = g*a(:,j+1); %formula for projection
    for k = 1:j
        Gn1(j+1) = Gn1(j+1) + F(k-1)*EQ*a(:,j-k+1);
    end
    noise1(j+1) = F(j)*b1;
    Gn1(j+1) = Gn1(j+1) + noise1(j+1);
end


%% approximated formula with limited history
a(:,1) = P*c;
b1 = Q*c;
for i = 1:n
    a(:,i+1) = EP*a(:,i);
    for k = 2:min(i,m)
        a(:,i+1) = a(:,i+1) + EP*(EQ^(k-1))*a(:,i-k+1);
    end
end

for j = 0:n
    Mn(j+1) = g*a(:,j+1); %formula for projection
    for k = 1:min(j,m)
        Mn(j+1) = Mn(j+1) + F(k-1)*EQ*a(:,j-k+1);
    end
end

%% Direct method
odefun = @(t,x) A*x;
[t,y]=ode113(odefun, tt ,c); % solve the equation
G = g*y(end,:)' % Gn = g*cn
G0 = g*y';

gg = E^n*g'; % gn = E^n*g
G1 = gg'*c % Gn = gn*c0

%% analytical solution
Gexact = g*(c.*exp(diag(A)*TT))

plot(tt,Gn,'b-o');
hold on;
plot(tt,Mn,'r-o');
hold on
plot(tt,G0,'k-o'); %I don't know why there is no lines through the point
