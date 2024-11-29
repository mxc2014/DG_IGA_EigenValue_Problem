function ref_nurbs = generate_ref_nurbs(nurbs_original, n_subdomains, Refinement)

ref_nurbs = cell(n_subdomains,1);

t = 3;

t_subdomains       = t*ones(n_subdomains,1);
Refine_subdomains  = Refinement*ones(n_subdomains,1);

e = 7; Refine_subdomains(e) = Refinement + 2;
e = 9; Refine_subdomains(e) = Refinement + 2;

for e = 1:n_subdomains
  
    [knotU,knotV] =  updateKnotVectors(nurbs_original{e}.knotU,nurbs_original{e}.knotV,nurbs_original{e}.a,nurbs_original{e}.b);%%  We want to have a regular mesh
    [knotU,knotV]=IGADegreeElevSurface(knotU,knotV,t_subdomains(e)); % Now we should consider the ratio of lengths in x- and y-directions
    pu = nurbs_original{e}.pu + t_subdomains(e);     pv = nurbs_original{e}.pv + t_subdomains(e); 
    
    ref_nurbs{e} = IGA_2D_Grid(knotU,knotV,pu,pv, Refine_subdomains(e));
    
    
end


for e = 2:n_subdomains % Update the the dofs information that is across the boundary of two sub-domains
    ref_nurbs{e} = Update_Edge_DoFs( ref_nurbs{e-1}.n_dofs_domains,ref_nurbs{e}   );
    
end


end