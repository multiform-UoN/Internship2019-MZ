noE2
evolprojfun = chebfun(@(y) evolprojected(odefun,P1,y,dt),[-1,1])
plot(evolortfun)
evolortfun = chebfun(@(y) evolorthogonal(odefun,P1,y,dt,1),[-1,1])
hold on
plot(evolortfun)
evolortfun2 = chebfun(@(y) evolorthogonal(odefun,P1,y,dt,2),[-1,1])