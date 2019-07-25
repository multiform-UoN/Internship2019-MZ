clear all
%close all

%This is the special case for g = P1
%% parameters
N = 100; % number of dimension
J = 1; % J < N
A = diag(1:N);
%A = rand(N);
A = -A; % coefficents for ODE
c = ones(N,1); % initial conditions
%% time
TT = 100;   % total time
dt = 1; % time step
tt = 0:dt:TT;  % time vector
n = length(tt)-1; % number of time steps
m = 6; % m < n number of past values to keep
E = expm(A*dt); % evolution matrix

tic

%% coeffcient
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

P = P1(:,1)*P1(:,1)'; 
for i = 2:J
    P = P + P1(:,i)*P1(:,i)';
end
P = P';
Q = eye(N) - P; % complementary matrix, orthogonal complement

%% precomputed evolution matrices
EQ = Q*E;
EP = P*E;
E1 = P1'*E*P1;
for k = 1:n
    E2(:,:,k) = P1'*E*EQ^k*P1;
end

%% Initialisation and time Loop
cpexact = zeros(J,n+1);
cp1 = zeros(J,n+1); % MZ filter
cp2 = zeros(J,n+1);
cpexact1 = zeros(J,n+1);% Truncated filter
noise1 = zeros(N,n+1);
cpexact(:,1) = P1'*c;
cp1(:,1) = P1'*c;
cp2(:,1) = P1'*c;
cpexact1(:,1) = P1'*c;
b1 = Q*c;

for i = 1:n
    cpexact(:,i+1) = E1*cpexact(:,i);
    for k = 2:i
        cpexact(:,i+1) = cpexact(:,i+1) + E2(:,:,k-1)*cpexact(:,i-k+1);
    end
    noise(:,i+1) = P1'*E*EQ^i*b1;
    cpexact(:,i+1) = cpexact(:,i+1) + noise(:,i+1);
end


for i = 1:n
    cp1(:,i+1) = E1*cp1(:,i);
    for k = 2:i
        cp1(:,i+1) = cp1(:,i+1) + E2(:,:,k-1)*cp1(:,i-k+1);
    end
end

for i = 1:n
    cp2(:,i+1) = E1*cp2(:,i);
    for k = 2:min(i,m)
        cp2(:,i+1) = cp2(:,i+1) + E2(:,:,k-1)*cp2(:,i-k+1);
    end
end

%% For J = 1
E11 = E1/(sum(E2)+E1);
for k = 1:m
    E21(:,:,k) = E2(:,:,k)/(sum(E2)+E1);
end
for i = 1:n
    cpexact1(:,i+1) = E11*cpexact1(:,i);
    for k = 2:min(i,m)
        cpexact1(:,i+1) = cpexact1(:,i+1) + E21(:,:,k-1)*cpexact1(:,i-k+1);
    end
end

figure;
error11 = abs(cpexact - cp1)%./abs(cpexact);
error21 = abs(cpexact - cp2)%./abs(cpexact);
error31 = abs(cpexact - cpexact1)%./abs(cpexact);
loglog(tt,error11,'-o');
hold on
plot(tt, error21,'-o');
plot(tt, error31,'-o');
xlabel('t');
ylabel('error');
legend('errors for m = n', 'errors for m = 6', 'Truncated filter');

figure;
error1 = abs(cpexact - cp1)./abs(cpexact);
error2 = abs(cpexact - cp2)./abs(cpexact);
error3 = abs(cpexact - cpexact1)./abs(cpexact);
loglog(tt,error1,'-o');
hold on
plot(tt, error2,'-o');
plot(tt, error3,'-o');
xlabel('t');
ylabel('error');
legend('errors for m = n', 'errors for m = 6', 'Truncated filter');


