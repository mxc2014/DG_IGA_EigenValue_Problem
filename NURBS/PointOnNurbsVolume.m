function S=PointOnNurbsVolume(ConPts,weights,knotU,pu,u,knotV,pv,v,knotW,pw,w)


ndim=size(ConPts,4);

Pw=WightedConPtsVolume(ConPts,weights);
Sw=PointOnBspVolume(Pw,knotU,pu,u,knotV,pv,v,knotW,pw,w);


S=Project(Sw);


end
 
