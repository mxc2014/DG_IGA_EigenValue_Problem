function DG_err   = Compute_DG_Error_Left_Right(nurbs_original_1,nurbs_refine_1,uh_1,nurbs_refine_2,uh_2, ref_nurbs_1, ref_uh_1,ref_nurbs_2, ref_uh_2)

% Now the two sub-domains are the left one  and the right one, the
% interface is at u = 1 for the left one, and u = 0 for the right one

% From the left sub-domain to the right sub-domain

% Here, the index 1 represents the left one, the index 2 represents the right one 
%%

ConPts_o_1   =  nurbs_original_1.ConPts;
weights_o_1  =  nurbs_original_1.weights;
knotU_o_1    =  nurbs_original_1.knotU; 
knotV_o_1    =  nurbs_original_1.knotV; 
pu_o_1       =  nurbs_original_1.pu;
pv_o_1       =  nurbs_original_1.pv;



knotV_1 = nurbs_refine_1.Vbar;
pu_1    = nurbs_refine_1.pu;
pv_1    = nurbs_refine_1.pv;


knotV_2 = nurbs_refine_2.Vbar;
pv_2    = nurbs_refine_2.pv;




%%




% right_edge_node_1 = nurbs_refine_1.right_edge_node;




right_edge_dofs_1 = nurbs_refine_1.right_edge_dofs_1st;
left_edge_dofs_2  = nurbs_refine_2.left_edge_dofs_1st;



VBreaks_1 = nurbs_refine_1.VBreaks;
VBreaks_2 = nurbs_refine_2.VBreaks;




VBreaks_merge = [VBreaks_1, VBreaks_2];
VBreaks_merge = unique(VBreaks_merge);
vNoEs = length(VBreaks_merge) - 1; 
 

right_edge_node_1 = zeros(vNoEs,2);

for i=1:vNoEs
right_edge_node_1(i,:) = [VBreaks_merge(i),VBreaks_merge(i+1)];
end




edge_element_1 =  Merge_edge_elements_idx(VBreaks_1,VBreaks_merge);
edge_element_2 =  Merge_edge_elements_idx(VBreaks_2,VBreaks_merge);


[gp,gw] = grule(pu_1+3); % The Gaussian quadrature rule on [-1,1].

n_gp = length(gp);


% As the interface is u = 1, so we consider u =1

% From the left sub-domain to the right sub-domain

left_interface_u  = 1; % The parameter u at the interface for the left sub-domain is 1
right_interface_u = 0; % The parameter u at the interface for the right sub-domain is 0


DG_err = 0;


     
for e=1:vNoEs
    
    
    local_DG_err = 0;
    
    edge_idx_plus  = edge_element_1(e);
    edge_idx_minus = edge_element_2(e);
    
    edge_dofs_plus  = right_edge_dofs_1(edge_idx_plus,:);
    edge_dofs_minus = left_edge_dofs_2(edge_idx_minus,:);
    

    
    ve = right_edge_node_1(e,:);
    a  = ve(1);    b = ve(2);
    J1 = (b-a)/2;  % The Jacobian from [-1,1] to [a,b].

    
    for i=1:n_gp
    v  = ((b-a)*gp(i) +a+b)/2;
    [~,DF_plus]   =  NurbsSurface(ConPts_o_1 ,weights_o_1,knotU_o_1 ,pu_o_1,left_interface_u, knotV_o_1,pv_o_1 ,v);% Left patch
    
    tau    = DF_plus(:,2);
    ds     = norm(tau);
   
       
    
    Jacobi = J1*gw(i)*ds; % The interface is u = 1, now along the v-direction
    
      
    Nu_plus     = 1;  % The non-zero basis fucntions along the interface u = 1 for the left sub-domain    
    Nv_plus     = bsplinebasis(knotV_1,pv_1,v);   % A column vector      
    basis_plus  = Nu_plus*Nv_plus;   
    uh_plus     = uh_1(edge_dofs_plus)'*basis_plus;
    ref_uh_plus = reference_uh(ref_nurbs_1,ref_uh_1,left_interface_u,v); 

    
        
    Nu_minus     = 1;  % The non-zero basis fucntions along the interface u = 0 for the right sub-domain    
    Nv_minus     = bsplinebasis(knotV_2,pv_2,v);   % A column vector      
    basis_minus  = Nu_minus*Nv_minus;  
    uh_minus     = uh_2(edge_dofs_minus)'*basis_minus;
    ref_uh_minus = reference_uh(ref_nurbs_2,ref_uh_2,right_interface_u,v); 
    
    local_DG_err = local_DG_err + (   (ref_uh_plus - uh_plus) - ( ref_uh_minus - uh_minus   )    ).^2*Jacobi;


   
    
   
    end
    
    
    DG_err   = DG_err  + local_DG_err;
    
    
end



 

end