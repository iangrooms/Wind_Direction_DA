function RHS = jacobian(q_hat,p)
% Function takes Fourier coefficients of vorticity (q_hat) and struct 
% containing parameters (p) and evaluates -fft[u dot grad vorticity]
% Jacobian is dealiased via 3/2 rule.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent DX DY
if isempty(DX)
    k = (2*pi/p.LX)*[0:.75*p.N-1 0 -.75*p.N+1:-1]';
    DX = 1i*repmat(k',[1.5*p.N 1]);
    DY = 1i*repmat(k,[1 1.5*p.N]);
    clear k
end

psi_hat = get_psi(q_hat,p);

% The following code computes the Jacobian J[psi,q], dealiased using a 3/2-rule.
% I.e., if you set N=256 then the code uses 3*N/2=384 Fourier modes (in each
% direction) to compute the Jacobian.
% factors of (4/9) and (9/4) are related to how Matlab normalizes the FFT
% on different size grids
Psi_hat = zeros(1.5*p.N);
Psi_hat(1:p.N/2+1,1:p.N/2+1) = (9/4)*psi_hat(1:p.N/2+1,1:p.N/2+1);
Psi_hat(1:p.N/2+1,p.N+2:1.5*p.N) = (9/4)*psi_hat(1:p.N/2+1,p.N/2+2:p.N);
Psi_hat(p.N+2:1.5*p.N,1:p.N/2+1) = (9/4)*psi_hat(p.N/2+2:p.N,1:p.N/2+1);
Psi_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N) = (9/4)*psi_hat(p.N/2+2:p.N,p.N/2+2:p.N);
Q_hat = zeros(1.5*p.N);
Q_hat(1:p.N/2+1,1:p.N/2+1) = (9/4)*q_hat(1:p.N/2+1,1:p.N/2+1);
Q_hat(1:p.N/2+1,p.N+2:1.5*p.N) = (9/4)*q_hat(1:p.N/2+1,p.N/2+2:p.N);
Q_hat(p.N+2:1.5*p.N,1:p.N/2+1) = (9/4)*q_hat(p.N/2+2:p.N,1:p.N/2+1);
Q_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N) = (9/4)*q_hat(p.N/2+2:p.N,p.N/2+2:p.N);
% calculate u.gradq on 3/2 grid
u = real(ifft2(-DY.*Psi_hat));
v = real(ifft2( DX.*Psi_hat));
qx= real(ifft2( DX.*Q_hat));
qy= real(ifft2( DY.*Q_hat));
jaco_real = u.*qx+v.*qy;
% fft, 3/2 grid; factor of (4/9) scales fft
Jaco_hat = (4/9)*fft2(jaco_real);
% reduce to normal grid
jaco_hat = zeros(p.N);
jaco_hat(1:p.N/2+1,1:p.N/2+1) = Jaco_hat(1:p.N/2+1,1:p.N/2+1);
jaco_hat(1:p.N/2+1,p.N/2+2:p.N) = Jaco_hat(1:p.N/2+1,p.N+2:1.5*p.N);
jaco_hat(p.N/2+2:p.N,1:p.N/2+1) = Jaco_hat(p.N+2:1.5*p.N,1:p.N/2+1);
jaco_hat(p.N/2+2:p.N,p.N/2+2:p.N) = Jaco_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N);

RHS = - jaco_hat;
