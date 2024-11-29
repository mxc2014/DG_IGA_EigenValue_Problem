function [M,rhs]= L2_Projection_Boundary_Front(nurbs_original, nurbs_refine,g, n_dofs)
%%
ConPts_o  = nurbs_original.ConPts;
weights_o = nurbs_original.weights;
knotU_o   = nurbs_original.knotU;
knotV_o   = nurbs_original.knotV;
knotW_o   = nurbs_original.knotW;
pu_o      = nurbs_original.pu;
pv_o      = nurbs_original.pv;
pw_o      = nurbs_original.pw;

knotV      = nurbs_refine.Vbar; 
knotW      = nurbs_refine.Wbar;  
pv         = nurbs_refine.pv;
pw         = nurbs_refine.pw; 
vNoEs      = nurbs_refine.vNoEs;
wNoEs      = nurbs_refine.wNoEs;

NoEs       = vNoEs * wNoEs; % The number of elements in the back face

VBreaks    = nurbs_refine.VBreaks; 
WBreaks    = nurbs_refine.WBreaks; 

coordinate = zeros(NoEs,4);

for i = 1:vNoEs
    for j = 1:wNoEs
        e = i + (j-1)*vNoEs;
        coordinate(e,:) = [VBreaks(i:(i+1)), WBreaks(j:(j+1))];
    end
end


ele_n_dofs = (pv+1)*(pw+1);

%%


front_face_dofs       = nurbs_refine.front_face_eles_1_layers_dofs ;


np      = pv + 1 ;
[gp,gw] = grule(np);

%%

u_front = 1.0;


M   = sparse(n_dofs,n_dofs);

rhs = zeros(n_dofs,1);


Me = zeros(ele_n_dofs,ele_n_dofs);
Fe = zeros(ele_n_dofs,1);

for e = 1:NoEs
    
    row  = front_face_dofs(e,:);
    Me = 0*Me; Fe = 0*Fe;
    ve = coordinate(e,1:2);  we = coordinate(e,3:4);
    
    va = ve(1);  vb = ve(2);  vJ = (vb - va)/2;
    wa = we(1);  wb = we(2);  wJ = (wb - wa)/2;
   
     

  for i = 1:np
  v = ( (vb - va )*gp(i) + va + vb  )/2;
  Nv = bsplinebasis(knotV,pv,v); % A column vector
  for j = 1:np
  w = ( (wb - wa )*gp(j) + wa + wb  )/2;
  Nw = bsplinebasis(knotW,pw,w); % A column vector
  
  N = Nv*Nw'; 
  N = N(:);

[F,DF] = NurbsVolume(ConPts_o,weights_o,knotU_o,pu_o,u_front,knotV_o,pv_o,v,knotW_o,pw_o,w);

A = DF(2,1)*DF(3,2) - DF(2,2)*DF(3,1);
B = DF(3,1)*DF(1,2) - DF(1,1)*DF(3,2);
C = DF(1,1)*DF(2,2) - DF(1,2)*DF(2,1);

J = sqrt(A*A+B*B+C*C)*vJ*wJ*gw(i)*gw(j);

Me = Me + N*N'*J;

Fe = Fe + g(F(1),F(2),F(3))*N*J;


end
   end
   


   M(row,row) = M(row,row) +  Me;
   rhs(row)   = rhs(row)   +  Fe;
   
end


end