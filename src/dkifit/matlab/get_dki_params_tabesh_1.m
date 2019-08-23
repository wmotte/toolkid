function dki_params=get_dki_params_tabesh_1(dvec,kvec,dki_tabs)

% Function that calculates the dki params MK,FK etc from kurtosis vector
% The vector is assumed to be of the form
% kvec is of the form k1111 k2222 k3333 k1112 k1113 k1222 k2223 k1333 k2333 k1122 k1133 k2233 k1123 k1223 k1233

temp=dvec([1 2 4 3 5 6]);
D_mtx=zeros(3);
D_mtx(1,1)=temp(1);
D_mtx(2,2)=temp(2);
D_mtx(3,3)=temp(3);
D_mtx(1,2)=temp(4)/2;% We assume that off-diagonal dvecs are actually a factor of 2 of the needed (to compensate for repeated terms)
D_mtx(2,1)=D_mtx(1,2);
D_mtx(2,3)=temp(5)/2;
D_mtx(3,2)=D_mtx(2,3);
D_mtx(1,3)=temp(6)/2;
D_mtx(3,1)=D_mtx(1,3);

D_mtx = [   1   0 0.5
            0   1 0
            0.5 0 1  ];

[V,D]=eig(D_mtx); % eigenvalues and eigenmatrix [3 x 3] from diffusion tensor

evals=diag(D); % trace
[ evals,indices ] = sort( evals,'descend' );
V = V(:, indices );

kvec_expanded=kvec(dki_tabs.kvec_to_table);

wtilda1111=V(dki_tabs.indices_wt1111(:,1)).*V(dki_tabs.indices_wt1111(:,2)).*V(dki_tabs.indices_wt1111(:,3)).*V(dki_tabs.indices_wt1111(:,4));
wtilda1111=transpose(wtilda1111(:))*kvec_expanded;

wtilda2222=V(dki_tabs.indices_wt2222(:,1)).*V(dki_tabs.indices_wt2222(:,2)).*V(dki_tabs.indices_wt2222(:,3)).*V(dki_tabs.indices_wt2222(:,4));
wtilda2222=transpose(wtilda2222(:))*kvec_expanded;

wtilda3333=V(dki_tabs.indices_wt3333(:,1)).*V(dki_tabs.indices_wt3333(:,2)).*V(dki_tabs.indices_wt3333(:,3)).*V(dki_tabs.indices_wt3333(:,4));
wtilda3333=transpose(wtilda3333(:))*kvec_expanded;

wtilda2233=V(dki_tabs.indices_wt2233(:,1)).*V(dki_tabs.indices_wt2233(:,2)).*V(dki_tabs.indices_wt2233(:,3)).*V(dki_tabs.indices_wt2233(:,4));
wtilda2233=transpose(wtilda2233(:))*kvec_expanded;

wtilda1133=V(dki_tabs.indices_wt1133(:,1)).*V(dki_tabs.indices_wt1133(:,2)).*V(dki_tabs.indices_wt1133(:,3)).*V(dki_tabs.indices_wt1133(:,4));
wtilda1133=transpose(wtilda1133(:))*kvec_expanded;

wtilda1122=V(dki_tabs.indices_wt1122(:,1)).*V(dki_tabs.indices_wt1122(:,2)).*V(dki_tabs.indices_wt1122(:,3)).*V(dki_tabs.indices_wt1122(:,4));
wtilda1122=transpose(wtilda1122(:))*kvec_expanded;

% MK
mk=f1_tabesh(evals)*wtilda1111+f1_tabesh([evals(2) evals(1) evals(3)])*wtilda2222+f1_tabesh([evals(3) evals(2) evals(1)])*wtilda3333+...
    f2_tabesh(evals)*wtilda2233+f2_tabesh([evals(2) evals(1) evals(3)])*wtilda1133+f2_tabesh([evals(3) evals(2) evals(1)])*wtilda1122;

%kpar
kpar=(sum(evals)^2/(9*(evals(1)^2)))*wtilda1111;

%kperp
kperp=g1_tabesh(evals)*wtilda2222+g1_tabesh([evals(1) evals(3) evals(2)])*wtilda3333+g2_tabesh(evals)*wtilda2233;

dki_params=[mk;kpar;kperp];

return;