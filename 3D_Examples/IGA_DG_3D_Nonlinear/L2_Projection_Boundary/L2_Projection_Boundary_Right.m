function [M,rhs]= L2_Projection_Boundary_Right(nurbs_original, nurbs_refine,g, n_dofs)
%%
ConPts_o  = nurbs_original.ConPts;
weights_o = nurbs_original.weights;
knotU_o   = nurbs_original.knotU;
knotV_o   = nurbs_original.knotV;
knotW_o   = nurbs_original.knotW;
pu_o      = nurbs_original.pu;
pv_o      = nurbs_original.pv;
pw_o      = nurbs_original.pw;

knotU      = nurbs_refine.Ubar; 
knotW      = nurbs_refine.Wbar;  
pu         = nurbs_refine.pu;
pw         = nurbs_refine.pw; 
uNoEs      = nurbs_refine.uNoEs;
wNoEs      = nurbs_refine.wNoEs;

NoEs       = uNoEs * wNoEs; % The number of elements in the left face

UBreaks    = nurbs_refine.UBreaks; 
WBreaks    = nurbs_refine.WBreaks; 
coordinate = zeros(NoEs,4);

for i = 1:uNoEs
    for j = 1:wNoEs
        e = i + (j-1)*uNoEs;
        coordinate(e,:) = [UBreaks(i:(i+1)), WBreaks(j:(j+1))];
    end
end


ele_n_dofs = (pu+1)*(pw+1);

%%

 

right_face_dofs       = nurbs_refine.right_face_eles_1_layers_dofs ;



np      = pu + 1 ;
[gp,gw] = grule(np);

%%

v_right = 1.0;


M   = sparse(n_dofs,n_dofs);

rhs = zeros(n_dofs,1);


Me = zeros(ele_n_dofs,ele_n_dofs);
Fe = zeros(ele_n_dofs,1);

for e = 1:NoEs
    
    row  = right_face_dofs(e,:);
    Me = 0*Me; Fe = 0*Fe;
    ue = coordinate(e,1:2);  we = coordinate(e,3:4);
    
    ua = ue(1);  ub = ue(2);  uJ = (ub - ua)/2;
    wa = we(1);  wb = we(2);  wJ = (wb - wa)/2;
   
     

  for i = 1:np
  u = ( (ub - ua )*gp(i) + ua + ub  )/2;
  Nu = bsplinebasis(knotU,pu,u); % A column vector
  for j = 1:np
  w = ( (wb - wa )*gp(j) + wa + wb  )/2;
  Nw = bsplinebasis(knotW,pw,w); % A column vector
  
  N = Nu*Nw'; 
  N = N(:);

[F,DF] = NurbsVolume(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v_right,knotW_o,pw_o,w);

A = DF(2,1)*DF(3,2) - DF(2,2)*DF(3,1);
B = DF(3,1)*DF(1,2) - DF(1,1)*DF(3,2);
C = DF(1,1)*DF(2,2) - DF(1,2)*DF(2,1);

J = sqrt(A*A+B*B+C*C)*uJ*wJ*gw(i)*gw(j);

Me = Me + N*N'*J;

Fe = Fe + g(F(1),F(2),F(3))*N*J;


end
   end
   


   M(row,row) = M(row,row) +  Me;
   rhs(row)   = rhs(row)   +  Fe;
   
end


end