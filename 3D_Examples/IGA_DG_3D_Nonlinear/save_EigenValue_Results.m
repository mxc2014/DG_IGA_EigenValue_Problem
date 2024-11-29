function save_EigenValue_Results(Lambda_h, Energy, n_eles, DOFs, p )

n_refine = size(Lambda_h,1);

fid = fopen( strcat('./EigenValue_Results/', 'Lambda_h_p_',num2str(p),'.txt'  ),'w');


fprintf( fid,'%s   p = %i   \n','# The fisrt four eigenvalues',p);

fprintf( fid,'No. of elements     || No. of DOFs    ||     The first eigenvalue   || Total energy\n');

for i = 1:n_refine
    fprintf(fid,'%i   %i   %1.15f   %1.15f \n', n_eles(i), DOFs(i), Lambda_h(i), Energy(i) );
end


end