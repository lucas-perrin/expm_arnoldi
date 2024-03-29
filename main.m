size_vec = 10;
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
fprintf(['convergence at default toler at n = ',num2str(n_bk),'\n'])
tic
expAb_matlab  = expm(A)*b;
time_matlab = toc;

error = norm(expAb_anroldi - expAb_matlab)/norm(expAb_anroldi);

fprintf('\n')
fprintf('size \\ times | Arnoldi  |  Matlab  |    Error   |  Speedup \n')
fprintf(['    ',num2str(prod(dimensions),'%6d'),'     | ',num2str(time_arnoldi,'%.6f'),' | ',num2str(time_matlab,'%.6f'),' | ',num2str(error,'%.4e'),' | ',num2str(time_matlab/time_arnoldi*100,'%.4f'),'%% \n\n'])