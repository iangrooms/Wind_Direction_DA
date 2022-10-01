Nt = 1024; % Number of assimilation cycles
load ../Reference/qp_ref_01.mat
ind_o = 16:32:256;
n_o = length(ind_o);
k = [0:128 -127:-1];
[KX,KY] = meshgrid(k,k);
K2 = KX.^2+KY.^2;
iK2 = 1./K2;iK2(1,1) = 0;
uK = 1i*KY./K2; uK(1,1) = 0;
vK =-1i*KX./K2; vK(1,1) = 0;
theta_obs = zeros([n_o n_o Nt]);
psi_obs = zeros([n_o n_o Nt]);
psi_obs_err = 0.002;
rng(1); % For reproducibility
for ii=1:Nt
    q_hat = fft2(qp_ref(:,:,ii));
    u = real(ifft2(uK.*q_hat));
    v = real(ifft2(vK.*q_hat));
    theta_obs(:,:,ii) = get_wind_direction(u(ind_o,ind_o),v(ind_o,ind_o));
    p = real(ifft2(-iK2.*q_hat));
    psi_obs(:,:,ii) = p(ind_o,ind_o) + psi_obs_err*randn(n_o);
end
clear u v uK vK p ind_x ind_y
