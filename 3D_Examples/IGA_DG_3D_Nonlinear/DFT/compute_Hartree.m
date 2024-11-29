function Hartree_Mat = compute_Hartree(nurbs_original,nurbs_refine,Hartree_h)

% Compute the the non-linear term

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

np = pu + 1;

[gp,gw] = grule(np);

%%







n_ele_dofs = (pu+1)*(pv+1)*(pw+1);

 
Hartree_e = zeros(n_ele_dofs,n_ele_dofs);


%  Max_nz = 1e9;
%  A = spalloc(n_dofs,n_dofs, Max_nz );






 Hartree_value  = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
 
 
 row_index = ones(NoEs*n_ele_dofs*n_ele_dofs,1);
 
 column_index = row_index;

 global_index = 1;

% tic


% bar = waitbar(0,'正在组装刚度矩阵');

for k=1:wNoEs
   wa = WBreaks(k); wb = WBreaks(k+1);  wJ = (wb-wa)/2; 
    for j=1:vNoEs
   va = VBreaks(j); vb = VBreaks(j+1);  vJ = (vb-va)/2; 
        for i=1:uNoEs
   ua = UBreaks(i); ub = UBreaks(i+1);  uJ = (ub-ua)/2;     

            Hartree_e =0*Hartree_e; 
            e = i + (j-1)*uNoEs + (k-1)*uNoEs*vNoEs;
            row = Element(e,:);

    for k1=1:np
    
    w     = ( (wb-wa)*gp(k1)+ wa + wb )/2;
     
    w_basis = bsplinebasis(knotW,pw,w); % A column vector
    
  
    for j1=1:np
      v     = ( (vb-va)*gp(j1)+ va + vb )/2;
      
      v_basis = bsplinebasis(knotV,pv,v); % A column vector
      
    for i1=1:np
      u     = ( (ub-ua)*gp(i1)+ ua + ub )/2;
      
      u_basis  = bsplinebasis(knotU,pu,u); % A column vector
      
      [~,DF]= NurbsVolume(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v,knotW_o,pw_o,w);
      
      vw_basis = v_basis*w_basis';  % A matrix
      vw_basis = vw_basis(:)';      % A row vector
      uvw_basis= u_basis*vw_basis;
      uvw_basis= uvw_basis(:);
    
      Hartree_in_ele = Hartree_h(row)'*uvw_basis;

      J = uJ*vJ*wJ*gw(k1)*gw(j1)*gw(i1)*abs(det(DF));
            
      Hartree_e = Hartree_e +  Hartree_in_ele*(uvw_basis*uvw_basis')*J;
     
 
      
end
end
end



             for i1=1:n_ele_dofs
                for j1=1:n_ele_dofs
                     row_index(global_index)    = row(i1);
                     column_index(global_index) = row(j1);
                      
                     Hartree_value(global_index) =  Hartree_e(i1,j1);
                     global_index = global_index + 1;
                end
             end
           
            
%             A(row,row) = A(row,row) + Ae;

% A = A + sparse(row(:),row(:),Ae(:),n_dofs,n_dofs);


           
                                           

        end
    end
end





 
  Hartree_Mat = sparse(row_index,column_index,Hartree_value,n_dofs,n_dofs); 
  
 clear row_index column_index value

 
%  whos row_index column_index  value  A
% 
%      % toc

end
