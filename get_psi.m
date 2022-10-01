function psi_hat = get_psi(q_hat,p)
persistent iDel
if(isempty(iDel))
    k = (2*pi/p.LX)*[0:p.N/2 -p.N/2+1:-1]';
    dX = 1i*repmat(k',[p.N 1]);
    dY = 1i*repmat(k,[1 p.N]);
    Laplacian = dX.^2+dY.^2;
    iDel = 1./Laplacian; iDel(1,1) = 0;
end
% Invert for psi
psi_hat = iDel.*q_hat;
