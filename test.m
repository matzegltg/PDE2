%%
gd = gd_3;
ns = ns_3;
sf = sf_3;
dl = decsg(gd,sf,ns);
[u,p,e,t] = adaptmesh(dl,@circleb1,1,0,0,"Maxt",500,"Ngen",Inf);
[p,e,t] = refinemesh(dl,p,e,t); 
pdemesh(p,e,t)