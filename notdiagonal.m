A = [1,2,3,4,5,6,0,8;1,2,4,3,0,6,7,8;1,0,3,4,5,6,8,7;2,1,0,4,5,6,7,8;0,2,3,4,5,7,6,8;1,3,2,0,5,6,7,8;1,2,3,0,6,5,7,8;3,2,1,4,5,0,7,8];
A = -A;
odefun = @(t,x) A*x;
c=ones(8,1);
tt=0:0.1:10;
[t,y]=ode113(odefun, tt ,c);
