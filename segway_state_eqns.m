function dx = segway_state_eqns(x, k1, k2, k3, k4, k5, k6, tpc, u, const_theta, holding, linear)

    if holding == true
        Sy = sin(const_theta);
        Cy = cos(const_theta);
        x(4) = 0;
        if linear == true
            Sy = const_theta;
            Cy = 1;
        end
    else
        Sy = sin(x(3));
        Cy = cos(x(3));
        if linear == true
            Sy = x(3);
            Cy = 1;
        end
    end

    dx(1,1) = x(2);
    dx(2,1) = (1/(k1*k4 - k2*k3*Cy^2))*((k4 + tpc*k2*Cy)*u + k6*k4*Sy*x(4)^2 - k5*k2*Sy*Cy);
    dx(3,1) = x(4);
    if holding == true
        dx(4,1) = 0;
    else
    dx(4,1) = (1/(k2*k3*Cy^2 - k1*k4))*((k3*Cy + tpc*k1)*u + k6*k3*Cy*Sy*x(4)^2 - k5*k1*Sy);
    end
end