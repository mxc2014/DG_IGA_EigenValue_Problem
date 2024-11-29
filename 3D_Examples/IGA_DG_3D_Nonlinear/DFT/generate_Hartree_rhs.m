function rhs = generate_Hartree_rhs(nurbs_original,nurbs_refine_domains, u0_h, rho, n_dofs,mp)

    rhs = zeros(n_dofs,1);
    
    n_subdomains = mp.n_subdomains;
    
    row = 1:nurbs_refine_domains{1}.n_dofs_domains;
    rhs(row) = generate_Hartree_rhs_subdomain(nurbs_original{1},nurbs_refine_domains{1}, u0_h(row), rho);
    
    for e = 2:n_subdomains
        row = ( 1 + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
        rhs(row) = generate_Hartree_rhs_subdomain(nurbs_original{e},nurbs_refine_domains{e}, u0_h(row), rho);     
    end
    
    
end