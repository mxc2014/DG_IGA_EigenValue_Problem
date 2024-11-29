function L2_norm = compute_L2_uh(nurbs_original,nurbs_refine,uh)

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;


Element=nurbs_refine.Element;
knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
 
Coordinate = nurbs_refine.Coordinate; 


Fhat = @(x,a,b) ( (b-a)*x+a+b )/2;
 

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;

 
 

pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;




n_ele_dofs = (pu+1)*(pv+1);


 
 



 
 


   L2_norm = 0;

 

       n_gp = pu + 3;
       
       [gp,gw] = grule(n_gp);

        for e = 1:NoEs
            
                     
            ue = Coordinate(e,1:2);   J1 = (ue(2) - ue(1) )/2;
            ve = Coordinate(e,3:4);   J2 = (ve(2) - ve(1) )/2;
            
 
            
            
                
            row = Element(e,:);
            
            for i=1:n_gp
                u       = Fhat(gp(i),ue(1),ue(2));
                Nu      = bsplinebasis(knotU,pu,u); % A column vector
               
                
                for j=1:n_gp
                v       = Fhat(gp(j),ve(1),ve(2));
                Nv      = bsplinebasis(knotV,pv,v);
                Nv      = Nv'; % A row vector
               
               
                
                [~,DF] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v);
                
 
                
                basis_funcs = Nu*Nv;  basis_funcs = basis_funcs(:);
                
                uh_at_gs_pnt = uh(row)'*basis_funcs;
                
                Jacobian = J1*J2*abs(det(DF))*gw(i)*gw(j);
                
 
                L2_norm = L2_norm + uh_at_gs_pnt.^2*Jacobian;   
                
%                  Fe = Fe + f(F(1,k1),F(2,k1))*uv_basis_ij(:,k1)*Jacobian(k1);
                end
            end
            

            

%            
            

%          n_gp = pu+3;                                  

        end
        
%         disp('r_min = ')
%         disp(r_min)

    

      
%       L2_norm = sqrt(L2_norm);
  

% disp('******************************************')
% disp('r_min =')
% disp(r_min)
% disp('******************************************')


end