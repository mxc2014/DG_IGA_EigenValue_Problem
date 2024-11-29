function DG_err   = Compute_DG_Error_Top_Bottom(nurbs_original_1,nurbs_refine_1,uh_1,nurbs_refine_2,uh_2, ref_nurbs_1, ref_uh_1,ref_nurbs_2, ref_uh_2)

% Now the two sub-domains are the upper one (v=0)  and lower one (v=1), the
% interface is at v = 0 for upper one, and v = 1 for the upper one

% From the bottom sub-domain to the upper sub-domain

% Here, the index 1 represents the bottom one, the index 2 represents the upper one 
%%

ConPts_o_1   =  nurbs_original_1.ConPts;
weights_o_1  =  nurbs_original_1.weights;
knotU_o_1    =  nurbs_original_1.knotU; 
knotV_o_1    =  nurbs_original_1.knotV; 
pu_o_1       =  nurbs_original_1.pu;
pv_o_1       =  nurbs_original_1.pv;


knotU_1 = nurbs_refine_1.Ubar;
pu_1    = nurbs_refine_1.pu;
 

knotU_2 = nurbs_refine_2.Ubar;
pu_2    = nurbs_refine_2.pu;
 
%%



% The plus sign points from the lower domain to upper element

UBreaks_1 = nurbs_refine_1.UBreaks;
UBreaks_2 = nurbs_refine_2.UBreaks;


top_edge_dofs_1     = nurbs_refine_1.top_edge_dofs_1st; % The dofs lying on the top edge of the lower subdomain
bottom_edge_dofs_2  = nurbs_refine_2.bottom_edge_dofs_1st;




UBreaks_merge = [UBreaks_1, UBreaks_2];
UBreaks_merge = unique(UBreaks_merge);
uNoEs         = length(UBreaks_merge) - 1; 



bottom_edge_node_1 = zeros(uNoEs,2);

for i=1:uNoEs
bottom_edge_node_1(i,:) = [UBreaks_merge(i),UBreaks_merge(i+1)];
end




edge_element_1 =  Merge_edge_elements_idx(UBreaks_1,UBreaks_merge);
edge_element_2 =  Merge_edge_elements_idx(UBreaks_2,UBreaks_merge);

[gp,gw] = grule(pu_1+3); % The Gaussian quadrature rule on [-1,1].

n_gp = length(gp);



v_top_1 = 1;     % The parameter of the top face of the bottom sub-domain

v_bottom_2 = 0;  % The parameter of the bottom face of the top sub-domain


% From lower sub-domain to upper sub-domain. It points from the lower to
% the upper.
% The plus is at the lower sub-domain, and minus is at the upper sub-domain


DG_err = 0;
     
for e=1:uNoEs
    
    edge_idx_plus  = edge_element_1(e);
    edge_idx_minus = edge_element_2(e);
    
    edge_dofs_plus  = top_edge_dofs_1(edge_idx_plus,:);
    edge_dofs_minus = bottom_edge_dofs_2(edge_idx_minus,:);
    
    
    
    ue = bottom_edge_node_1(e,:);
    a  = ue(1);    b = ue(2);
    J1 = (b-a)/2;  % The Jacobian from [-1,1] to [a,b].
  
    
    local_DG_err = 0;
    
    for i=1:n_gp
    u  = ((b-a)*gp(i) +a+b)/2;
    [~,DF_plus]   =  NurbsSurface(ConPts_o_1 ,weights_o_1,knotU_o_1 ,pu_o_1,u,knotV_o_1,pv_o_1,v_top_1);% Left patch
    
    tau    = DF_plus(:,1);
    ds     = norm(tau);
     
    
   
    Jacobi = J1*gw(i)*ds; % The interface is u = 1, now along the v-direction
    
    Nu_plus     =  bsplinebasis(knotU_1,pu_1,u);  % A column vector  
    Nv_plus     =  1;        
    basis_plus  = Nu_plus*Nv_plus;   % basis_plus = basis_plus(:);
    uh_plus     = uh_1(edge_dofs_plus)'*basis_plus;
    ref_uh_plus = reference_uh(ref_nurbs_1,ref_uh_1,u,v_top_1); 
    

        
   
    Nu_minus     = bsplinebasis(knotU_2,pu_2,u);   % For Omega 2, it is the bottom boundary, v = 0.   
    Nv_minus     = 1;               
    basis_minus  = Nu_minus*Nv_minus;   % basis_minus = basis_minus(:);
    uh_minus     = uh_2(edge_dofs_minus)'*basis_minus;
    ref_uh_minus = reference_uh(ref_nurbs_2,ref_uh_2,u,v_bottom_2);

    
    local_DG_err = local_DG_err + ( ref_uh_plus -  uh_plus  - (  ref_uh_minus - uh_minus  )  ).^2*Jacobi;
    
   
   
    end
    
    
    DG_err = DG_err + local_DG_err;
    
    
     
   
    
end



 

end