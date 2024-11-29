function M_rho  = compute_Nonlinear_Mass_Mat(nurbs_original,nurbs_refine,uh,rho)

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


 
 Me =zeros(n_ele_dofs,n_ele_dofs); % Element mass matrix



% Fe = zeros(n_ele_dofs,1);

% A = sparse(n_dofs,n_dofs);
% rhs = zeros(n_dofs,1);



 
M_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
 

row_index    = M_value;
column_index = M_value;
global_index = 1;





%        A = sparse(n_dofs,n_dofs); 
%        M =sparse(n_dofs,n_dofs); 
%        M_Vext = sparse(n_dofs,n_dofs);


       
       
%         r_min = 100;

       n_gp = pu + 3;
       
       [gp,gw] = grule(n_gp);

        for e = 1:NoEs
            
            Me = 0*Me;               
            ue = Coordinate(e,1:2);   J1 = (ue(2) - ue(1) )/2;
            ve = Coordinate(e,3:4);   J2 = (ve(2) - ve(1) )/2;
            
%             if (ue(1) - 0.5)*(ue(2) - 0.5)<0 && (ve(1) - 0.5)*(ve(2) - 0.5)<0
%                   disp(ue)
%                   disp(ve)
%                 if mod(pu,2)==1  % pu is an odd integer
%                    n_gp = pu + 100;
%                 else
%                    n_gp = pu + 100;
%                 end
%             end
            
%             [gp,gw] = grule(n_gp);
            
            
            
                
            row = Element(e,:);
            
            for i=1:n_gp
                u       = Fhat(gp(i),ue(1),ue(2));
                Nu      = bsplinebasis(knotU,pu,u); % A column vector
               
                
                for j=1:n_gp
                v       = Fhat(gp(j),ve(1),ve(2));
                Nv      = bsplinebasis(knotV,pv,v);
                Nv      = Nv'; % A row vector
               
               
                
                [~,DF] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v);
                
%                 r = norm(F);
                
%                 r_min = min(r,r_min); 
%                 if r<=1.0e-5
%                      disp('the quadrature point is near singular point')
%                      disp(F)
%                      disp(n_gp)      
%                      r = 1.0e-3;
%                 end
                
                basis_funcs = Nu*Nv;  basis_funcs = basis_funcs(:);
                
                uh_at_gs_pnt = uh(row)'*basis_funcs;
                
                Jacobian = J1*J2*abs(det(DF))*gw(i)*gw(j);
                
 
                Me = Me + rho(uh_at_gs_pnt)*(basis_funcs*basis_funcs')*Jacobian;     
%                  Fe = Fe + f(F(1,k1),F(2,k1))*uv_basis_ij(:,k1)*Jacobian(k1);
                end
            end
            

            
             for i1=1:n_ele_dofs
                for j1=1:n_ele_dofs
                     row_index(global_index) = row(i1);
                     column_index(global_index) = row(j1);
                     M_value(global_index) =  Me(i1,j1);       
                     global_index = global_index + 1;

%                     A(row,row) = A(row,row) + Ae;
%                     M(row,row) = M(row,row) + Me;
%                     M_Vext(row,row) = M_Vext(row,row)  +  Me_ext;
                    
                 end
             end
%            
            

%          n_gp = pu+3;                                  

        end
        
%         disp('r_min = ')
%         disp(r_min)

    

      
       M_rho   =  sparse(row_index,column_index,M_value,n_dofs,n_dofs); 
  

% disp('******************************************')
% disp('r_min =')
% disp(r_min)
% disp('******************************************')


end