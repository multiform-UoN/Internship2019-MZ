function Fp = evolF(odefun,g,P1,y,dt,i)
y1 = P1 * y;
Fp1 = y1;
for j=1:i
    Fp1 = evol(odefun,Fp1,dt)';
    Fp1 = Fp1 - P1*P1'*Fp1;
end
Fp = g * Fp1;
end
