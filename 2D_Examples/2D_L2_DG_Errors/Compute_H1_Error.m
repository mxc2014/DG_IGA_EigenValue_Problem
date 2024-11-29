function err   = Compute_H1_Error(nurbs_original,nurbs_refine,uh,ref_uh,ref_nurbs, ref_nurbs_original)

% This script is used to compute the $L^2$ norm, $H^1$ semi-norm and $H^2$
% semi-norm  error between $u$ and the IGA solution $u_h$.



ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;


pu   =  nurbs_refine.pu;
pv   =  nurbs_refine.pv;




Element=nurbs_refine.Element;
Coordinate = nurbs_refine.Coordinate;


Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
 

 

NoEs=nurbs_refine.NoEs;







 u_np = pu + 1 ;   % It seems that if pu >=3, u_np = pu is enough!!!
 v_np = pv + 1 ;   
 
% u_np = ceil((pu+1)/2) + 1;
% 
% v_np = ceil((pu+1)/2) + 1;


[gp_u,gw_u]=grule(u_np); % The quadrature points and quadrature weights in the [-1,1].

gp_v = gp_u;
gw_v = gw_u;



L2_err = 0;  
H1_err = 0;



for e=1:NoEs % Loop for each element;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);

    U_ij = uh(row);
    U_ij = U_ij';

    
    for j=1:v_np
            vJ = (ve(2) - ve(1))/2;
            v = vJ*gp_v(j) + (ve(1)+ve(2))/2;
            Vders=bspbasisDers(Vbar,pv,v,1);
            Nv=Vders(1,:); DNv=Vders(2,:);  
        for i=1:u_np
            uJ = (ue(2) - ue(1))/2;
            u = uJ*gp_u(i)+(ue(1)+ue(2))/2;
            Uders=bspbasisDers(Ubar,pu,u,1);
            Nu=Uders(1,:)';   DNu=Uders(2,:)';  
            
            [~,~, ~,~,DF]=NurbsSurfaceDers(ConPts_o,knotU_o,knotV_o,weights_o,pu_o,u,pv_o,v);
            
            [ref_uh_at_uv,grad_ref_uh_at_uv]   = grad_reference_uh(ref_nurbs, ref_nurbs_original, ref_uh,u,v);
            

            
            J=abs(det(DF))*gw_u(i)*uJ*gw_v(j)*vJ;
           
            B= Nu*Nv;
            B= B(:);
            
            L2_err = L2_err + ( ref_uh_at_uv - U_ij*B ).^2*J;
            
            
           DBu=DNu*Nv; DBu=DBu(:);
           DBv=Nu*DNv; DBv=DBv(:);

           DB=[DBu,DBv]/DF;


           H1_err = H1_err + sum(( grad_ref_uh_at_uv - U_ij*DB ).^2)*J;
               

        end
    end
    
    
            
end


%     L2_err  = sqrt(L2_err );
%     H1_err  = sqrt(H1_err );
   
    
    err = zeros(1,2);
    
    err(1) = L2_err;
    err(2) = H1_err;
     
end