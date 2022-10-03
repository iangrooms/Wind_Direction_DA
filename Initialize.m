N = 256; % Number of points in each direction
LX = 2*pi; % Nondimensional domain width
rDrag = 1E-2; % Coefficient of linear drag/friction
nu2 = 4E-4; % Coefficient of Laplacian vorticity diffusion
nu4 = 0*5E-8; % Coefficient of hyperviscous vorticity diffusion
beta = 10;

% Set up linear dissipation operator & forcing spectrum
k = (2*pi/LX)*[0:N/2 -N/2+1:-1]'; % wavenumbers
L = zeros(N);
S = L; % stochastic forcing, not spectrum though
for jj=1:N
    for ii=1:N
        kr = sqrt(k(ii)^2+k(jj)^2);
        L(ii,jj) = -(rDrag - beta*1i*k(jj)/(kr^2)+ nu2*kr^2 + nu4*kr^4);
        if( abs(kr - 11) <= 1)
            S(ii,jj) = 7.97e3;
        end
    end
end
L(1,1) = 0;
clear k kr ii jj

% Put useful stuff into a struct
params = struct('N',N,'LX',LX,'S',S,'L',L);

% Initialize
t = 0;
% qp is the "physical" vorticity, i.e. values on the grid
% Don't initialize to zero
qp = 1E-4*randn(N);
% q is the Fourier coefficients of vorticity
q = fft2(qp);
