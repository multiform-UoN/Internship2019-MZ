function E1 = evolprojected(odefun,P1,y,dt)
y1 = P1 * y;
Ey1 = evol(odefun,y1,dt)';
E1 = P1' * Ey1;
end

