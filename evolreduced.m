function EP1 = evolreduced(odefun,P1,dt)
tt = 0:dt:dt;
[a,b] = size(P1);
for i = 1:b
    EP1(:,i) = evol(odefun,P1(:,i),dt);
end
end