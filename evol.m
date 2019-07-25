function EZ = evol(odefun,z,dt)
tt = 0:dt:dt;
[t,y]=ode45(odefun, tt, z);
EZ = y(end, :);
end