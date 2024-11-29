function e_total = get_total_energy(nurbs_original,nurbs_refine,A,uh,rho,V_ext,Hartree_h,mp,uh_idx)

occ = 2;
e_kinetic = occ*uh'*A*uh/2;
e_ext     = 0;
e_hartree = 0;
e_xc      = 0;
e_total   = e_kinetic;

n_subdomains = mp.n_subdomains;

for e = 1:n_subdomains
    row = uh_idx{e};
    [e_total_subDomain,e_ext_subDomain,e_hartree_subDomain,e_xc_subDomain] =  get_energy_subDomain(nurbs_original{e},nurbs_refine{e},uh(row),rho,V_ext,Hartree_h(row));
    e_total   = e_total   + e_total_subDomain;
    e_ext     = e_ext     + e_ext_subDomain;
    e_hartree = e_hartree + e_hartree_subDomain;
    e_xc      = e_xc      + e_xc_subDomain;
    
end

% disp('**************************************************************')
% disp('e_total e_ext e_hartree e_xc')
% disp([e_total, e_ext, e_hartree, e_xc])
% disp('**************************************************************')
 







end