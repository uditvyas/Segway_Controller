function dx = segway_state_eqns(x, k1, k2, k3, k4, k5, u, const_theta, holding)
    
    if holding == true
        Sy = sin(const_theta);
        Cy = cos(const_theta);
        x(4) = 0;
    else
        Sy = sin(x(3));
        Cy = cos(x(3));
    end
    dx(1,1) = x(2);
    dx(2,1) = (1/(k1*k4 - k2*k3*Cy))*((k4 + k2*Cy)*u + k2*k4*Sy*x(4)^2 - k5*k2*Sy*Cy);
    dx(3,1) = x(4);
    dx(4,1) = (1/(k2*k3*Cy - k1*k4))*((k3*Cy + k1)*u + k2*k3*Cy*Sy*x(4)^2 - k5*k1*Sy);
end