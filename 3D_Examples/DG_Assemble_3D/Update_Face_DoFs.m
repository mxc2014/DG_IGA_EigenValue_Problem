function nurbs_refine = Update_Face_DoFs(n_dofs_1,nurbs_refine)

% From Omega_1 to Omega_2, we have to update the index of dofs on Omega_2

% nurbs_refine.n_dofs = n_dofs_1 +  nurbs_refine.n_dofs;

nurbs_refine.bottom_face_dof                =  n_dofs_1 +  nurbs_refine.bottom_face_dof ;
nurbs_refine.bottom_face_eles_2layers_dofs  =  n_dofs_1 +  nurbs_refine.bottom_face_eles_2layers_dofs  ;
nurbs_refine.bottom_face_eles_1_layers_dofs =  n_dofs_1 +  nurbs_refine.bottom_face_eles_1_layers_dofs  ;


nurbs_refine.top_face_dof                   =  n_dofs_1 +  nurbs_refine.top_face_dof;
nurbs_refine.top_face_eles_2layers_dofs     =  n_dofs_1 +  nurbs_refine.top_face_eles_2layers_dofs   ;
nurbs_refine.top_face_eles_1_layers_dofs    =  n_dofs_1 +  nurbs_refine.top_face_eles_1_layers_dofs   ;


nurbs_refine.front_face_dof                   =  n_dofs_1 +  nurbs_refine.front_face_dof;
nurbs_refine.front_face_eles_2layers_dofs     =  n_dofs_1 +  nurbs_refine.front_face_eles_2layers_dofs   ;
nurbs_refine.front_face_eles_1_layers_dofs    =  n_dofs_1 +  nurbs_refine.front_face_eles_1_layers_dofs   ;

nurbs_refine.back_face_dof                    =  n_dofs_1  +   nurbs_refine.back_face_dof; 
nurbs_refine.back_face_eles_2layers_dofs      =  n_dofs_1  +   nurbs_refine.back_face_eles_2layers_dofs    ;
nurbs_refine.back_face_eles_1_layers_dofs     =  n_dofs_1  +   nurbs_refine.back_face_eles_1_layers_dofs    ;

nurbs_refine.left_face_dof                    =  n_dofs_1  +  nurbs_refine.left_face_dof;
nurbs_refine.left_face_eles_2layers_dofs      =  n_dofs_1  +  nurbs_refine.left_face_eles_2layers_dofs     ;
nurbs_refine.left_face_eles_1_layers_dofs     =  n_dofs_1  +  nurbs_refine.left_face_eles_1_layers_dofs     ;

nurbs_refine.right_face_dof                   =  n_dofs_1  +  nurbs_refine.right_face_dof; 
nurbs_refine.right_face_eles_2layers_dofs     =  n_dofs_1  +  nurbs_refine.right_face_eles_2layers_dofs    ;
nurbs_refine.right_face_eles_1_layers_dofs    =  n_dofs_1  +  nurbs_refine.right_face_eles_1_layers_dofs    ;


end