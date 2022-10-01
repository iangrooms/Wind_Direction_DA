% Runs a series of DA cycles using categorical wind direction observations

% Allocate storage for forecast/analysis mean and spread
FM = zeros([256 256 Nt]);
FS = FM;
AM = FM;
AS = FS;

% Allocate storage for wind ensemble and fft of streamfunction
p = zeros(Ne,1);
psi_hat = 1i*zeros([256 256 Ne]);

LR = zeros(8,8,Nt);

% Run the cycled DA experiment
for ii=1:Nt
    % Forecast mean and spread
    FM(:,:,ii) = mean(Q,3);
    FS(:,:,ii) = std(Q,0,3);
    fprintf('Forecast RMSE %1.2f, ',sqrt(mean(mean((FM(:,:,ii)-qp_ref(:,:,ii)).^2))))
    obs_mask = zeros(256);
    for ki=1:n_o
        for kj=1:n_o
            obs_mask(ind_o(ki), ind_o(kj)) = 1;
        % Assimilate p
            % Get obs ensemble
            for n=1:Ne
                psi_hat(:,:,n) = get_psi(fft2(Q(:,:,n)),params);
                tmp = real(ifft2(psi_hat(:,:,n)));
                p(n) = tmp(ind_o(ki),ind_o(kj));
            end
            % Ensemble perturbation matrix
            A = reshape(Q-mean(Q,3),[],Ne)/sqrt(Ne-1);
            % Ensemble obs perturbation matrix, p
            V = (p - mean(p))'/sqrt(Ne-1);
            s2 = V*(V'); g2 = psi_obs_err^2; % prior and likelihood variances in obs space
            % Set up localization
            p_inc = (s2/(s2+g2))*psi_obs(ki,kj,ii) + (psi_obs_err^2/(s2+g2))*mean(p) ...
                    + (psi_obs_err/sqrt(s2+g2))*(p - mean(p));
            p_inc = reshape(p_inc - p,[1 1 Ne]);
            rho = max(abs(p_inc(:)));
            locRad1 = locRad*(rho+.0064)/(rho+0.0256);
            LR(ki,kj,ii) = locRad1;
            loc = exp(-(.5/locRad1^2)*(x.^2 + x'.^2));
            loc1 = circshift(loc,[ind_o(ki)-1 ind_o(kj)-1]);
            % Update mean
            Q = Q + reshape(bsxfun(@times,loc1(:),(A*(V'))*(psi_obs(ki,kj,ii)-mean(p))*(1/(s2 + g2))),[256 256]);
            % Update perturbations
            WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2)); % Whitaker & Hamil Beta
            Q = Q - WHB*sqrt(Ne-1)*reshape(bsxfun(@times,loc1(:),A*(V')*V),[256 256 Ne]);
        end
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
    parfor(jj=1:Ne, 18)
        Q(:,:,jj) = real(ifft2(forecast(fft2(Q(:,:,jj)),params)));
    end
    if( mod(ii,32)==0 )
        fprintf('Finished cycle %04d \n',ii)
    end
end
