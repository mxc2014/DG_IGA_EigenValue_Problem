function Hartree_Mat = compute_Hartree_XC(nurbs_original,nurbs_refine,Hartree_h, rho, uh)

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

 
Hartree_XC_e = zeros(n_ele_dofs,n_ele_dofs);


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

            Hartree_XC_e =0*Hartree_XC_e; 
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
      
      rho_uh_e = rho( uh(row)'* uvw_basis );
      
      %%
      
      rs      = ((3/(4*pi))./rho_uh_e).^(1/3);
      falpha = -0.458165293283143;
      vx = 4/3*falpha./rs;
      ux = falpha./rs; 
      
      
      
a = 0.0311;
b = -0.048;
c = 0.0020;
d = -0.0116;
gc = -0.1423;
b1 = 1.0529;
b2 = 0.3334;

% vc = zeros(size(rs));
% uc = zeros(size(rs));
% 
% % high density formula
% idxl = rs < 1;
% rsl = rs(idxl);
% lnrs = log (rsl);
% uc(idxl) = a*lnrs + b + c*rsl.*lnrs + d*rsl;
% vc(idxl) = a*lnrs + (b-a/3) + 2/3*c*rsl.*lnrs + (2*d-c)/3*rsl;
% % interpolation formula
% idxg = rs >= 1;
% rsg = rs(idxg);
% rs12 = sqrt(rsg);
% ox = 1 + b1*rs12 + b2*rsg;
% dox = 1 + 7/6*b1*rs12 + 4/3*b2*rsg;
% uc(idxg) = gc./ox;
% vc(idxg) = uc(idxg).*dox./ox;


if rs < 1
      lnrs = log (rs);
      uc = a*lnrs + b + c*rs*lnrs + d*rs;
      vc = a*lnrs + (b-a/3) + 2/3*c*rs*lnrs + (2*d-c)/3*rs;
    
else
% interpolation formula
      rs12 = sqrt(rs);
      ox   = 1 + b1*rs12 + b2*rs;
      dox  = 1 + 7/6*b1*rs12 + 4/3*b2*rs;
      uc   = gc./ox;
      vc   = uc.*dox./ox;
end

vxc = vx+vc;
uxc = ux+uc;



      
      %%
            
      Hartree_XC_e = Hartree_XC_e +  (Hartree_in_ele + vxc )*(uvw_basis*uvw_basis')*J;
     
 
      
end
end
end



             for i1=1:n_ele_dofs
                for j1=1:n_ele_dofs
                     row_index(global_index)    = row(i1);
                     column_index(global_index) = row(j1);
                      
                     Hartree_value(global_index) =  Hartree_XC_e(i1,j1);
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
