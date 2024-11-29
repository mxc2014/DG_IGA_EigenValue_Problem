function [A, V_ext, M] = solve_laplace_A_M(nurbs_original,nurbs_refine,V_ext)

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
knotW_o   = nurbs_original.knotW; 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;
pw_o        =  nurbs_original.pw;




Element=nurbs_refine.Element;
Coordinate = nurbs_refine.Coordinate;





knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
knotW=nurbs_refine.Wbar;
 

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;



pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;
pw            =  nurbs_refine.pw;



%% Gauss quadrature rule

np = pu + 1;
[gp,gw] = grule(np);

%%







n_ele_dofs = (pu+1)*(pv+1)*(pw+1);

Ae        = zeros(n_ele_dofs,n_ele_dofs);
Me        = zeros(n_ele_dofs,n_ele_dofs);
V_ext_e   = zeros(n_ele_dofs,n_ele_dofs);









 A_value    = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
 Vext_value = A_value;
 M_value    = A_value;
 
 row_index = ones(NoEs*n_ele_dofs*n_ele_dofs,1);
 column_index = row_index;
 global_index = 1;

% tic


% bar = waitbar(0,'正在组装刚度矩阵');

for e = 1:NoEs
           ua = Coordinate(e,1); ub = Coordinate(e,2);  uJ = (ub-ua)/2;
           va = Coordinate(e,3); vb = Coordinate(e,4);  vJ = (vb-va)/2;
           wa = Coordinate(e,5); wb = Coordinate(e,6);  wJ = (wb-wa)/2;   

            Ae = 0*Ae;  V_ext_e = 0*V_ext_e; Me =0*Me; 
            
            row = Element(e,:);

    for k1=1:np
    w       = ( (wb-wa)*gp(k1)+ wa + wb )/2;
    wDers   = bspbasisDers(knotW,pw,w,1);
    w_basis = wDers(1,:);
    Dw_basis = wDers(2,:);
    for j1=1:np
      v     = ( (vb-va)*gp(j1)+ va + vb )/2;
      vDers = bspbasisDers(knotV,pv,v,1);
      v_basis = vDers(1,:);
      Dv_basis = vDers(2,:);
    for i1=1:np
      u     = ( (ub-ua)*gp(i1)+ ua + ub )/2;
      uDers = bspbasisDers(knotU,pu,u,1);
      u_basis  = uDers(1,:);
      Du_basis = uDers(2,:);
      [F,DF]= NurbsVolume(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v,knotW_o,pw_o,w);
      
      vw_basis = v_basis'*w_basis; vw_basis = vw_basis(:)';
      uvw_basis= u_basis'*vw_basis;
      uvw_basis= uvw_basis(:);
    
      Du_vw_basis = Du_basis'*vw_basis;   Du_vw_basis = Du_vw_basis(:);
 
      Dv_w_basis  = Dv_basis'*w_basis; Dv_w_basis = Dv_w_basis(:)';
      Dv_uw_basis = u_basis'*Dv_w_basis;  Dv_uw_basis = Dv_uw_basis(:);
      
      Dw_v_basis   = v_basis'*Dw_basis; Dw_v_basis = Dw_v_basis(:)';
      Dw_uv_basis = u_basis'*Dw_v_basis;   Dw_uv_basis=Dw_uv_basis(:);
      
      basis_grad = [Du_vw_basis,Dv_uw_basis,Dw_uv_basis]/DF;

      J = uJ*vJ*wJ*gw(k1)*gw(j1)*gw(i1)*abs(det(DF));
      
      temp =  (uvw_basis*uvw_basis')*J;

      Ae      = Ae + basis_grad*basis_grad'*J ;
      V_ext_e = V_ext_e  + V_ext(F(1),F(2),F(3))*temp;
      Me      = Me + temp;
     
 
      
end
end
end



             for i1=1:n_ele_dofs
                for j1=1:n_ele_dofs
                     row_index(global_index)    = row(i1);
                     column_index(global_index) = row(j1);
                     A_value(global_index)    =   Ae(i1,j1);
                     Vext_value(global_index) =   V_ext_e(i1,j1);
                     M_value(global_index)    =   Me(i1,j1);
                     global_index = global_index + 1;
                end
             end
           
            
%             A(row,row) = A(row,row) + Ae;

% A = A + sparse(row(:),row(:),Ae(:),n_dofs,n_dofs);


end





  A     = sparse(row_index,column_index,A_value,n_dofs,n_dofs);    
  V_ext = sparse(row_index,column_index,Vext_value,n_dofs,n_dofs); 
  M     = sparse(row_index,column_index,M_value,n_dofs,n_dofs); 
  
 clear row_index column_index value

 
%  whos row_index column_index  value  A
% 
%      % toc

end
