% Runs a series of DA cycles using categorical wind direction observations

% Allocate storage for forecast/analysis mean and spread
FM = zeros([256 256 Nt]);
FS = FM;
AM = FM;
AS = FS;

% Allocate storage for wind ensemble and fft of streamfunction
u = zeros(Ne,1);
v = zeros(Ne,1);
p = zeros(Ne,1);
psi_hat = 1i*zeros([256 256 Ne]);

LR = zeros(64,Nt);

% Run the cycled DA experiment
for ii=1:Nt
    % Forecast mean and spread
    FM(:,:,ii) = mean(Q,3);
    FS(:,:,ii) = std(Q,0,3);
    fprintf('Forecast RMSE %1.2f, ',sqrt(mean(mean((FM(:,:,ii)-qp_ref(:,:,ii)).^2))))
    for ko=1:n_o
        ind_ox = ind_o(ko,1,ii);
        ind_oy = ind_o(ko,2,ii);
        % Step 1: Get p, u, v increments
        % Get prior ensemble
        for n=1:Ne
            psi_hat(:,:,n) = get_psi(fft2(Q(:,:,n)),params);
            tmp = real(ifft2(psi_hat(:,:,n)));
            p(n) = tmp(ind_ox,ind_oy);
            tmp = real(ifft2(-1i*KY.*psi_hat(:,:,n)));
            u(n) = tmp(ind_ox,ind_oy);
            tmp = real(ifft2( 1i*KX.*psi_hat(:,:,n)));
            v(n) = tmp(ind_ox,ind_oy);
        end
        % Get p increments
        s2 = var(p);g2 = psi_obs_err^2;
        p_inc = (s2/(s2+g2))*psi_obs(ko,ii) + (psi_obs_err^2/(s2+g2))*mean(p) ...
                + (psi_obs_err/sqrt(s2+g2))*(p - mean(p));
        p_inc = reshape(p_inc - p,[1 1 Ne]);
        rho = max(abs(p_inc(:)));
        locRad1 = locRad*(rho+.0064)/(rho+0.0256);
        LR(ko,ii) = locRad1;
        % Get u & v increments
        % Get angles relative to observed
        theta = mod(atan2(v,u) - (pi/16)*theta_obs(ko,ii),2*pi);
        theta(theta>pi) = theta(theta>pi) - 2*pi;
        % Get prior pdf estimate using KDE
        [prior_pdf,R] = get_prior_pdf(theta);
        % Get prior cdf
        prior_cdf = (pi/128)*cumtrapz(prior_pdf);
        prior_cdf(end+1) = 1;
        % Get posterior pdf
        posterior_pdf = l.*prior_pdf;
        posterior_pdf = posterior_pdf/( (pi/128)*sum(posterior_pdf) );
        % Get posterior cdf
        posterior_cdf = (pi/128)*cumtrapz(posterior_pdf(117:139));
        posterior_cdf(end+1) = 1;
        % Transform prior thetas to uniform
        theta_u = interp1(theta_grid,prior_cdf,theta);
        % Apply inverse of posterior cdf to theta_u
        theta = interp1(posterior_cdf,theta_grid(116:139),theta_u);
        % Shift back to 0 being east
        theta = theta + (pi/16)*theta_obs(ko,ii);
        % Get observation ensemble increments (p, u, v)
        r = sqrt(u.^2+v.^2);
        u_inc = reshape(r.*cos(theta) - u,[1 1 Ne]);
        v_inc = reshape(r.*sin(theta) - v,[1 1 Ne]);
    % Step 2: Conditional sampling via regression
        % Set up localization
        loc = exp(-(.5/locRad1^2)*(x.^2 + x'.^2));
        loc1 = circshift(loc,[ind_ox-1 ind_oy-1]);
        % Use p, u, v to get regression coefficients 
        Z = [ones(Ne,1) p(:) u(:) v(:)];
        beta = shiftdim(reshape(Z\reshape(shiftdim(Q,2),Ne,[]),[4 256 256]),1);
        % Find linear regression increments
        Q_inc = bsxfun(@times,beta(:,:,2),p_inc)...
              + bsxfun(@times,beta(:,:,3),u_inc)...
              + bsxfun(@times,beta(:,:,4),v_inc);
        Q = Q + bsxfun(@times,loc1,Q_inc);
    end
    % Store EnKF analysis mean
    AM(:,:,ii) = mean(Q,3);
    % RTPS 20CRv3 used.3 + .7* in the NH and .1 + .9* in the SH
    AS(:,:,ii) = std(Q,0,3);
    rInf = (1-RTPS) + RTPS*(FS(:,:,ii)./AS(:,:,ii));
    Q = AM(:,:,ii) + bsxfun(@times,rInf,Q - AM(:,:,ii));
    % Store EnKF posterior spread
    AS(:,:,ii) = std(Q,0,3);
    fprintf('Analysis RMSE %1.2f\n',sqrt(mean(mean((AM(:,:,ii)-qp_ref(:,:,ii)).^2))))
    if(any(isnan(Q(:))))
        break
    end
    % forecast ensemble
    parfor jj=1:Ne
        Q(:,:,jj) = real(ifft2(forecast(fft2(Q(:,:,jj)),params)));
    end
    if( mod(ii,32)==0 )
        fprintf('Finished cycle %04d \n',ii)
    end
end
