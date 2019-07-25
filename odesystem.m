syms c1(t) c2(t) c3(t) c4(t) c5(t) c6(t) c7(t) c8(t)
A = [1,0,0,0,0,0,0,0;0,2,0,0,0,0,0,0;0,0,3,0,0,0,0,0;0,0,0,4,0,0,0,0;0,0,0,0,5,0,0,0;0,0,0,0,0,6,0,0;0,0,0,0,0,0,7,0;0,0,0,0,0,0,0,8];
A = -A;
C = [c1;c2;c3;c4;c5;c6;c7;c8];
odes = diff(C) == 0.1*A*C;
[c1Sol(t),c2Sol(t),c3Sol(t),c4Sol(t),c5Sol(t),c6Sol(t),c7Sol(t),c8Sol(t)] = dsolve(odes);
c1Sol(t) = simplify(c1Sol(t));
c2Sol(t) = simplify(c2Sol(t));
c3Sol(t) = simplify(c3Sol(t));
c4Sol(t) = simplify(c4Sol(t));
c5Sol(t) = simplify(c5Sol(t));
c6Sol(t) = simplify(c6Sol(t));
c7Sol(t) = simplify(c7Sol(t));
c8Sol(t) = simplify(c8Sol(t));

C = C(0) == [1;1;1;1;1;1;1;1];
[c1Sol(t), c2Sol(t), c3Sol(t), c4Sol(t), c5Sol(t), c6Sol(t), c7Sol(t), c8Sol(t)] = dsolve(odes,C);

C1 = [c1Sol(t), c2Sol(t), c3Sol(t), c4Sol(t), c5Sol(t), c6Sol(t), c7Sol(t), c8Sol(t)];

i = 1:8;
for j = 0:10
c(j+1,:) = [c1Sol(j),c2Sol(j),c3Sol(j),c4Sol(j),c5Sol(j),c6Sol(j),c7Sol(j),c8Sol(j)];
end

for k = 1:11
plot(i,c(k,:),'-o');
hold on
end
grid on
legend('t=0','t=1','t=2','t=3','t=4','t=5','t=6','t=7','t=8','t=9','t=10');
xlabel('i');
ylabel('cis');