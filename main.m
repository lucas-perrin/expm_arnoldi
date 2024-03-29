size_vec = 15;
dim = 3;
dimensions = size_vec*ones(1,dim);

% 3D Dirichlet Boundary conditions laplacian on mesh of size [size_vec, size_vec, size_vec]
[lambda,V,A] = laplacian(dimensions); 

% random vector b
b = rand(prod(dimensions),1);

% compute exp(A)*b with Arnoldi and with expm and compare results
tic
[Q,h,expAb_anroldi,n_bk] = Arnoldi_exp(A,b,size_vec,1,40);
time_arnoldi = toc;
n_bk
tic
expAb_matlab  = expm(A)*b;
time_matlab = toc;

error = norm(expAb_anroldi - expAb_matlab)/norm(expAb_anroldi);

fprintf('\n')
fprintf('size \\ times | Arnoldi  |  Matlab   |   Error \n')
fprintf(['    ',num2str(prod(dimensions),'%6d'),'     | ',num2str(time_arnoldi,'%.6f'),' | ',num2str(time_matlab,'%.6f'),' | ',num2str(error,'%4e'),' \n\n'])