function [A,M] = solve_laplace_A_F(nurbs_original,nurbs_refine,V_ext_1)

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

 np = 2*pu + 6;

[gp,gw] = grule(np);

n_ele_dofs = (pu+1)*(pv+1);

Ae =zeros(n_ele_dofs,n_ele_dofs); % Element stiffness matrix
Me =zeros(n_ele_dofs,n_ele_dofs); % Element mass matrix
Me_ext =zeros(n_ele_dofs,n_ele_dofs); % Element mass matrix for V ext


Fe = zeros(n_ele_dofs,1);

% A = sparse(n_dofs,n_dofs);
% rhs = zeros(n_dofs,1);



A_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
M_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
% M_ext_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);

row_index = A_value;
column_index = A_value;
global_index = 1;





%        A = sparse(n_dofs,n_dofs); 
%        M =sparse(n_dofs,n_dofs); 
%        M_Vext = sparse(n_dofs,n_dofs);


       
       
%        r_min = 100;



        for e = 1:NoEs
            Ae = 0*Ae;  Me = 0*Me;  Me_ext = 0*Me_ext;    Fe =0*Fe;              
            ue = Coordinate(e,1:2);   J1 = (ue(2) - ue(1) )/2;
            ve = Coordinate(e,3:4);   J2 = (ve(2) - ve(1) )/2;
            row = Element(e,:);
            
            for i=1:np
                u       = Fhat(gp(i),ue(1),ue(2));
                uders   = bspbasisDers(knotU,pu,u,1);
                Nu  = uders(1,:)';
                DNu = uders(2,:)';
                
                for j=1:np
                v       = Fhat(gp(j),ve(1),ve(2));
                vders   = bspbasisDers(knotV,pv,v,1);
                Nv  = vders(1,:);
                DNv = vders(2,:);
                
                [F,DF] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v);
                
                basis_funcs = Nu*Nv;  basis_funcs = basis_funcs(:);
                DNu_v = DNu*Nv; DNu_v = DNu_v(:);
                DNv_u = Nu*DNv; DNv_u = DNv_u(:);
                
                basis_grad = [DNu_v,DNv_u]/DF;
                
                Jacobian = J1*J2*abs(det(DF))*gw(i)*gw(j);
                

                
                 
                
                
                Mass_mat_ele = basis_funcs*basis_funcs'*Jacobian; 
                
                Me = Me + Mass_mat_ele;
                
                Ae = Ae + basis_grad*basis_grad'*Jacobian/2 + ( V_ext_1(F(1),F(2)) )*Mass_mat_ele;
                
%                 r = norm(F(:,k1) - [1/4;1/2] );
                
%                 r_min = min(r,r_min);
                
                
%                 if(r<=1.0e-3)
%                     disp(r)
%                 end
                
%                 Me_ext  =  Me_ext + ( V_ext_1(F(1),F(2)) )*Mass_mat_ele;
                 
                end
            end
            

            
            for i1=1:n_ele_dofs
               for j1=1:n_ele_dofs
                    row_index(global_index) = row(i1);
                    column_index(global_index) = row(j1);
                    A_value(global_index) =  Ae(i1,j1);
                    M_value(global_index) =  Me(i1,j1);
%                     M_ext_value(global_index) =  Me_ext(i1,j1);                   
                    global_index = global_index + 1;


                end
            end
            
            




%                     A(row,row) = A(row,row) + Ae;
%                     M(row,row) = M(row,row) + Me;
%                     M_Vext(row,row) = M_Vext(row,row)  +  Me_ext;
                    
                                           

        end

    

      A =sparse(row_index,column_index,A_value,n_dofs,n_dofs); 
      M =sparse(row_index,column_index,M_value,n_dofs,n_dofs); 
%       M_Vext = sparse(row_index,column_index,M_ext_value,n_dofs,n_dofs);


% disp('******************************************')
% disp('r_min =')
% disp(r_min)
% disp('******************************************')


end