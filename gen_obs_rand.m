% Get obs
n_o = 64;
ind_o = zeros([n_o 2 Nt]); 
k = [0:128 -127:-1];
[KX,KY] = meshgrid(k,k);
K2 = KX.^2+KY.^2;
iK2 = 1./K2;iK2(1,1) = 0;
uK = 1i*KY./K2; uK(1,1) = 0;
vK =-1i*KX./K2; vK(1,1) = 0;
theta_obs = zeros([n_o Nt]);
psi_obs = zeros([n_o Nt]);
psi_obs_err = 0.002;
rng(1); % For reproducibility
for ii=1:Nt
    q_hat = fft2(qp_ref(:,:,ii));
    u = real(ifft2(uK.*q_hat));
    v = real(ifft2(vK.*q_hat));
    p = real(ifft2(-iK2.*q_hat));
    tmp = randperm(256);
    ind_o(:,1,ii) = tmp(1:n_o);
    tmp = randperm(256);
    ind_o(:,2,ii) = tmp(1:n_o);
    for jj=1:n_o
        theta_obs(jj,ii) = get_wind_direction(u(ind_o(jj,1,ii),ind_o(jj,2,ii)),...
                                              v(ind_o(jj,1,ii),ind_o(jj,2,ii)));
        psi_obs(jj,ii) = p(ind_o(jj,1,ii),ind_o(jj,2,ii)) + psi_obs_err*randn(1);
    end
end
clear u v uK vK p tmp
