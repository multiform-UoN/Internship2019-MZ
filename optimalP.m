A = rand(10);
A = A'*A + eye(10);
[V,D] = eig(A); % A = V*D*V'
P1 = V;
P = P1 * P1';


P11 = P1(:,1:4);
PP = P11*P11';
QQ = eye(10) - PP;
I1 = P11'*V;
I2 = V'*P11;



A = rand(10);
A = A'*A + eye(10);
[L,D] = ldl(A); % A = V*D*V'
P1 = L;
PP = P1 * P1';
QQ = eye(10) - PP;


P11 = P1(:,1:4);
PP = P11*P11';
QQ = eye(10) - PP;
I1 = P11'*V;
I2 = V'*P11;



