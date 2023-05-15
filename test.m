%%
gd = gd_3;
ns = ns_3;
sf = sf_3;
dl = decsg(gd,sf,ns);
[u,p,e,t] = adaptmesh(dl,"cirsb",1,0,0);

%%
c45 = cos(pi/4);
L1 = [2 -c45 0  c45 0 1 0 0 0 0]';
L2 = [2 -c45 0 -c45 0 1 0 0 0 0]';
C1 = [1 -c45  c45 -c45 -c45 1 0 0 0 1]';
C2 = [1  c45  c45 -c45  c45 1 0 0 0 1]';
C3 = [1  c45 -c45  c45  c45 1 0 0 0 1]';
g = [L1 L2 C1 C2 C3];

[u,p,e,t] = adaptmesh(g,"cirsb",1,0,0);
pdemesh(p,e,t)