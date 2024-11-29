function  [e_total,e_ext,e_hartree,e_xc] =  get_energy_subDomain(nurbs_original,nurbs_refine,uh,rho,V_ext,Hartree_h)



e_ext     = 0;
e_hartree = 0;
e_xc      = 0;

 



ConPts_o   =  nurbs_original.ConPts ;
weights_o  =  nurbs_original.weights;
knotU_o    =  nurbs_original.knotU; 
knotV_o    =  nurbs_original.knotV; 
knotW_o    =  nurbs_original.knotW; 
pu_o       =  nurbs_original.pu;
pv_o       =  nurbs_original.pv;
pw_o       =  nurbs_original.pw;




Element=nurbs_refine.Element;
Coordinate = nurbs_refine.Coordinate;





knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
knotW=nurbs_refine.Wbar;
 

NoEs=nurbs_refine.NoEs;
 



pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;
pw            =  nurbs_refine.pw;

%%

np = pu + 1;

[gp,gw] = grule(np);

%%

for e = 1:NoEs
           ua = Coordinate(e,1); ub = Coordinate(e,2);  uJ = (ub-ua)/2;
           va = Coordinate(e,3); vb = Coordinate(e,4);  vJ = (vb-va)/2;
           wa = Coordinate(e,5); wb = Coordinate(e,6);  wJ = (wb-wa)/2;   

           
            
            row = Element(e,:);

    for k1=1:np
    w          = ( (wb-wa)*gp(k1)+ wa + wb )/2;
    w_basis    = bsplinebasis(knotW,pw,w);
 
     
    for j1=1:np
      v       = ( (vb-va)*gp(j1)+ va + vb )/2;
      v_basis = bsplinebasis(knotV,pv,v);
     
      
    for i1=1:np
      u       = ( (ub-ua)*gp(i1)+ ua + ub )/2;
      u_basis = bsplinebasis(knotU,pu,u);
   
       
      [F,DF]= NurbsVolume(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v,knotW_o,pw_o,w);
      
      vw_basis = v_basis*w_basis';    vw_basis = vw_basis(:);
      uvw_basis= u_basis*vw_basis';   uvw_basis= uvw_basis(:);
    


      J = uJ*vJ*wJ*gw(k1)*gw(j1)*gw(i1)*abs(det(DF));
      
      uh_at_uvw    = uh(row)'*uvw_basis;
      rho_h_at_uvw = rho(uh_at_uvw);
      
      e_ext = e_ext +  V_ext(F(1),F(2),F(3))*rho_h_at_uvw*J;
      
      e_hartree = e_hartree + Hartree_h(row)'*uvw_basis*rho_h_at_uvw/2*J;
      
      
%%
      
      rs      = ((3/(4*pi))./rho_h_at_uvw).^(1/3);
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
e_xc = e_xc + uxc*rho_h_at_uvw*J;
      
 

      
%       V_ext_e = V_ext_e  + V_ext(F(1),F(2),F(3))*temp;
%       Me      = Me + temp;
     
 
      
end
end
    end
end

    e_nn = 0;
    
    e_total = e_ext + e_hartree + e_xc + e_nn;




end