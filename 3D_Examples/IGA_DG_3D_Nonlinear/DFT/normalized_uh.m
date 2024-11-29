function uh  = normalized_uh(nurbs_original, nurbs_refine_domains, uh, n_subdomains)

row = 1:nurbs_refine_domains{1}.n_dofs_domains;
L2_norm_uh = compute_L2_uh_subdomain( nurbs_original{1},nurbs_refine_domains{1}, uh(row) );

for e = 2:n_subdomains
   row = ( 1 + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
   L2_norm_uh =  L2_norm_uh + compute_L2_uh_subdomain(nurbs_original{e},nurbs_refine_domains{e}, uh(row) );
end

L2_norm_uh = sqrt(L2_norm_uh);
uh = uh/L2_norm_uh;


end