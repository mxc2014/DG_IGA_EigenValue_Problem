function [uh_at_uv,grad_uh_at_uv]   = grad_reference_uh(nurbs_refine,nurbs_original, uh,u,v)

% This script is used to compute the reference solution uh at the
% parametric point (u,v)


ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o   = nurbs_original.knotU; 
knotV_o   = nurbs_original.knotV; 
pu_o      =  nurbs_original.pu;
pv_o      =  nurbs_original.pv;


pu   =  nurbs_refine.pu;
pv   =  nurbs_refine.pv;


Ubar = nurbs_refine.Ubar;
Vbar = nurbs_refine.Vbar;

m = nurbs_refine.m;

i = findspan(Ubar,pu,u);
j = findspan(Vbar,pv,v);

n_ele_dofs = (pu+1)*(pv+1);

row = zeros(n_ele_dofs,1);
k = 1;

for j1 = (j-pv):j
    for i1 = (i-pu):i
         row(k) = i1 + (j1-1)*m;
         k = k+1;
    end
end


   Uders = bspbasisDers(Ubar,pu,u,1);
   Nu = Uders(1,:)'; DNu = Uders(2,:)'; % Now they are column vectors 
   
   Vders=bspbasisDers(Vbar,pv,v,1);
   Nv = Vders(1,:); DNv = Vders(2,:); 
   
   N = Nu*Nv;  N = N(:);
   
   uh_at_uv = uh(row)'*N;
   
   
   [~,~, ~,~,DF]=NurbsSurfaceDers(ConPts_o,knotU_o,knotV_o,weights_o,pu_o,u,pv_o,v);
   
   
    DBu=DNu*Nv; DBu=DBu(:);
    DBv=Nu*DNv; DBv=DBv(:);
    DB=[DBu,DBv]/DF;
    
    grad_uh_at_uv = uh(row)'*DB;
        


     
end