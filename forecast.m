function q_out = forecast(q,params)
% Integrates 2D Navier Stokes for nt time steps of size dt from initial
% condition q. q and q_out are Fourier coefficients.
L = params.L;
S = params.S;

t = 0;
nt = 0;
dt = 0.02;
tol = 1E-3;
r0 = tol;
tf = 0.1;

while t < tf
    M = 1./(1-.25*dt*L);
    % First stage ARK4
    k0 = jacobian(q,params);
    l0 = L.*q;
    % Second stage
    q1 = M.*(q+.5*dt*k0+.25*dt*l0);
    k1 = jacobian(q1,params);
    l1 = L.*q1;
    % Third stage
    q2 = M.*(q+dt*(13861*k0/62500+6889*k1/62500+8611*l0/62500-1743*l1/31250));
    k2 = jacobian(q2,params);
    l2 = L.*q2;
    % Fourth stage
    q3 = M.*(q+dt*(-0.04884659515311858*k0-0.1777206523264010*k1+0.8465672474795196*k2...
    +0.1446368660269822*l0-0.2239319076133447*l1+0.4492950415863626*l2));
    k3 = jacobian(q3,params);
    l3 = L.*q3;
    % Fifth stage
    q4 = M.*(q+dt*(-0.1554168584249155*k0-0.3567050098221991*k1+1.058725879868443*k2...
    +0.3033959883786719*k3+0.09825878328356477*l0-0.5915442428196704*l1...
    +0.8101210538282996*l2+0.2831644057078060*l3));
    k4 = jacobian(q4,params);
    l4 = L.*q4;
    % Sixth stage
    q5 = M.*(q+dt*(0.2014243506726763*k0+0.008742057842904184*k1+0.1599399570716811*k2...
    +0.4038290605220775*k3+0.2260645738906608*k4+0.1579162951616714*l0...
    +0.1867589405240008*l2+0.6805652953093346*l3-0.2752405309950067*l4));
    k5 = jacobian(q5,params);
    l5 = L.*q5;
    % Error control
    r1 = dt*max(max(max(abs(ifft2(0.003204494398459*(k0+l0) -0.002446251136679*(k2+l2)-0.021480075919587*(k3+l3)...
       +0.043946868068572*(k4+l4) -0.023225035410765*(k5+l5))))));
    if r1>tol
        % Reject step
        dt=.75*dt;
    else
        % Accept step
        nt = nt + 1;
        t = t + dt;
        qp = real(ifft2(q+dt*(0.1579162951616714*(k0+l0)+0.1867589405240008*(k2+l2)+...
        0.6805652953093346*(k3+l3)-0.2752405309950067*(k4+l4)+(k5+l5)/4)));
        q = fft2(qp);
        % Add some stochastic forcing (not during RK4 since not smooth)
        q = q + sqrt(dt)*S.*exp(1i*2*pi*rand(params.N));
        % step size adjustment: EPS, PI.3.4 ; divide by 4 for a 4th order
        % method with 3rd order embedded
        dt = min(((.75*tol/r1)^(.3/4))*((r0/r1)^(.4/4))*dt,tf-t);
        r0=r1;
    end
end
q_out = q;
%fprintf('Took %03d steps\n',nt)
