%Transfer matrix function for calculating reflection coefficients
function coeff = TransferMatrix(lambda, theta, n, d)

    k = (2*pi)/lambda;
    kz = n(1).*k.*sin(theta);
    kx = zeros(1,length(n));
    
    kx(1) = n(1)*k*cos(theta);
    for m = 2:length(n)
        kx(m) = sqrt((k*n(m))^2 - kz^2);
    end

    phi = zeros(1,length(n) - 2);
    for m = 1:length(d)
        phi(m) = kx(m + 1)*d(m);
    end
    
    M = eye(2);
    for m = 1:length(n) - 2
        A = 0.5*(1 + (kx(m + 1)./kx(m)));
        B =  0.5*(1 - (kx(m + 1)./kx(m)));
        D = [[A B];[B A]];
        P = [[exp(-1i*phi(m)) 0];[0 exp(1i*phi(m))]];
        M = M*D*P;
    end
    
    A = 0.5*(1 + (kx(end)./kx(end-1)));
    B =  0.5*(1 - (kx(end)./kx(end-1)));
    D = [[A B];[B A]];
    M = M*D;
    
    coeff = zeros(1,2);
    coeff(1) = M(2,1)/M(1,1);
    %coeff(1) = M(2,1)/M(1,1);
    coeff(2) = 1 - coeff(1);
    %coeff(2) = (abs(1/M(1,1)).^2)*real(k_3x/k_1x);
end
