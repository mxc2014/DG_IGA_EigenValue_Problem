function rhs = generate_Hartree_rhs_subdomain(nurbs_original,nurbs_refine, uh, rho)


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
UBreaks=nurbs_refine.UBreaks;   % u 方向上节点向量中的断�?
VBreaks=nurbs_refine.VBreaks;    % v 方向上节点向量中的断�?
WBreaks=nurbs_refine.WBreaks;    % w 方向上节点向量中的断�?


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

rhs = zeros(n_dofs,1);

Fe = zeros(n_ele_dofs,1);



for k=1:wNoEs
   wa = WBreaks(k); wb = WBreaks(k+1);  wJ = (wb-wa)/2; 
    for j=1:vNoEs
   va = VBreaks(j); vb = VBreaks(j+1);  vJ = (vb-va)/2; 
        for i=1:uNoEs
   ua = UBreaks(i); ub = UBreaks(i+1);  uJ = (ub-ua)/2;     

            Fe =0*Fe; 
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
      uvw_basis= uvw_basis(:);      % A column vector
    
      uh_in_ele = uh(row)'*uvw_basis;

      J = uJ*vJ*wJ*gw(k1)*gw(j1)*gw(i1)*abs(det(DF));
      Fe = Fe + 4*pi*rho(uh_in_ele)*uvw_basis*J;
          
end
end
    end

 rhs(row) = rhs(row) + Fe;


           
                                           

        end
    end
end


  
 clear row_index column_index value

 



end