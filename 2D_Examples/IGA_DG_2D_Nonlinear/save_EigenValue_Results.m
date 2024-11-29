function save_EigenValue_Results(Lambda_h,n_eles, DOFS, p )

n_refine = size(Lambda_h,1);

fid = fopen( strcat('./EigenValue_Results/', 'Lambda_h_p_',num2str(p),'.txt'  ),'w');


fprintf( fid,'%s   p= %i \n','# The fisrt eigenvalue',p );

fprintf( fid,'No. of elements     || No. of DOFs    ||     The first four eigenvalues\n\n');

for i = 1:n_refine
%     fprintf(fid,'%i   %i   %1.15f %1.15f %1.15f %1.15f\n', n_eles(i), DOFS(i), Lambda_h(i,1),Lambda_h(i,2),Lambda_h(i,3),Lambda_h(i,4) );
 fprintf(fid,'%i   %i   %1.15f  \n', n_eles(i), DOFS(i), Lambda_h(i,1) );
end


end