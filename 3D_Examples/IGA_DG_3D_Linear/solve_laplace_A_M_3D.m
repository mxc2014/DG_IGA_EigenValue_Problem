function [A,M_ext, M] = solve_laplace_A_M_3D(nurbs_original,nurbs_refine,V_ext,e)

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
knotW_o   = nurbs_original.knotW; 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;
pw_o        =  nurbs_original.pw;




Element=nurbs_refine.Element;





knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
knotW=nurbs_refine.Wbar;
UBreaks=nurbs_refine.UBreaks;   % u 方向上节点向量中的断点.
VBreaks=nurbs_refine.VBreaks;    % v 方向上节点向量中的断点.
WBreaks=nurbs_refine.WBreaks;    % w 方向上节点向量中的断点.

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;
wNoEs=nurbs_refine.wNoEs;

pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;
pw            =  nurbs_refine.pw;



%%

if e == 14
np = pu + 5;
else 
    np  = pu + 3;
end

[gp,gw] = grule(np);

%%







n_ele_dofs = (pu+1)*(pv+1)*(pw+1);

Ae = zeros(n_ele_dofs,n_ele_dofs);
Ve = Ae;
Me = zeros(n_ele_dofs,n_ele_dofs);

A = sparse(n_dofs,n_dofs);

M_ext = A;

M = sparse(n_dofs,n_dofs);


%  Max_nz = 1e9;
%  A = spalloc(n_dofs,n_dofs, Max_nz );






%  A_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
%  M_value = A_value;
%  row_index = ones(NoEs*n_ele_dofs*n_ele_dofs,1);
%  column_index = row_index;
%  global_index = 1;

% tic


% bar = waitbar(0,'正在组装刚度矩阵');

for k=1:wNoEs
   wa = WBreaks(k); wb = WBreaks(k+1);  wJ = (wb-wa)/2; 
    for j=1:vNoEs
   va = VBreaks(j); vb = VBreaks(j+1);  vJ = (vb-va)/2; 
        for i=1:uNoEs
   ua = UBreaks(i); ub = UBreaks(i+1);  uJ = (ub-ua)/2;     

            Ae = 0*Ae;  Ve = 0*Ve;  Me =0*Me; 
            e = i + (j-1)*uNoEs + (k-1)*uNoEs*vNoEs;
            row = Element(e,:);

    for k1=1:np
      w     = ( (wb-wa)*gp(k1)+ wa + wb )/2;
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

      Ae = Ae + basis_grad*basis_grad'*J/2 ;
      Ve = Ve + V_ext(F(1),F(2),F(3))*temp;
      
      Me = Me + temp;
     
 
      
end
end
end



%              for i1=1:n_ele_dofs
%                 for j1=1:n_ele_dofs
%                      row_index(global_index) = row(i1);
%                      column_index(global_index) = row(j1);
%                      A_value(global_index) = A_value(global_index) + Ae(i1,j1);
%                      M_value(global_index) = M_value(global_index) + Me(i1,j1);
%                      global_index = global_index + 1;
%                 end
%              end
           
            
            A(row,row) = A(row,row) + Ae;
            M_ext(row,row) = M_ext(row,row) + Ve;
            M(row,row) = M(row,row) + Me;

%  A = A + sparse(row(:),row(:),Ae(:),n_dofs,n_dofs);


           
                                           

        end
    end
end





%   A = sparse(row_index,column_index,A_value,n_dofs,n_dofs);    
%   A = A/2;
 
%   M = sparse(row_index,column_index,M_value,n_dofs,n_dofs); 
  
%  clear row_index column_index value

 
%  whos row_index column_index  value  A
% 
%      % toc

end
