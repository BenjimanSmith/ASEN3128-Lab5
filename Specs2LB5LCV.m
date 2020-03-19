% linear ODE Function
%Benjiman Smith
%02/13/2020
function dydt = Specs2LB5LCV(t, Conditions, Force, Disturb, givens)
    % define givens
    m = givens(6);      % mass of the drone [kg]
    r = givens(7);       % body to motor distance [m]
    k = givens(8);     % [Nm/N]
    R = givens(9);  % [m]
    g = givens(10);% gravity [m/s^2]
    alpha = givens(1);   % [N/(m/s)^2]
    eta = givens(2);    % [N/(rad/s)^2]
    I_x = givens(3);   % moment of inertia in the x direction [kg m^2]
    I_y = givens(4);   % moment of inertia in the y direction[kg m^2]
    I_z = givens(5);  % moment of inertia in the z direction[kg m^2]
    k_1 = givens(11);
    k_2 = givens(12);
    k_3 = givens(13);
    k_4 = givens(14);
    k_5 = givens(15);
    k_6 = givens(16);
        
    % define conditions vector
    deltau = Conditions(1); % inertial velocity in the u direction in body coordinates [m/s]
    deltav = Conditions(2); % inertial velocity in the v direction in body coordinates [m/s]
    deltaw = Conditions(3); % inertial velocity in the w direction in body coordinates [m/s]
    deltap = Conditions(4); % inertial angular velocity in the p direction in body coordinates [rad/s]
    deltaq = Conditions(5); % inertial angular velocity in the q direction in body coordinates [rad/s]
    deltar = Conditions(6); % inertial angular velocity in the r direction in body coordinates [rad/s]
   deltaphi = Conditions(7); % bank [rad]
    deltatheta = Conditions(8); % elevation [rad]
    psi = Conditions(9); % azimuth [rad]
    x_E = Conditions(10); % position vector in the x direction in inertial coordinates [m]
    y_E = Conditions(11); % position vector in the y direction in inertial coordinates [m]
    z_E = Conditions(12); % position vector in the z direction in inertial coordinates [m]
    f1 = Force(1); % trim force exerted by motor 1 [N]
    f2 = Force(2); % trim force exerted by motor 1 [N]
    f3 = Force(3); % trim force exerted by motor 1 [N]
    f4 = Force(4); % trim force exerted by motor 1 [N]
     if t <= 2
         deltavr = .5; % 2.1 m/s
         deltaur = 0;
     else
         deltavr = 0;
         deltaur=  0;
     end
    % Aerodynamic/Control Moments and Forces
    Laero = - alpha^2 * deltap^2 * sign(deltap); % p component of the aerodynamic moments
    Maero = - alpha^2 * deltaq^2 * sign(deltaq); % q component of the aerodynamic moments
    Naero = - alpha^2 * deltar^2 * sign(deltar); % w component of the aerodynamic moments
 Lcontrol = ((f1 + f2) - (f3 + f4)) * R; % p component of the control moments
    Mcontrol = ((f3 + f2) - (f1 + f4)) * R; % q component of the control moments
    Ncontrol = -0.004*deltar; % w component of the aerodynamic moments
    deltaL = Laero + Lcontrol; % L moment sum (abt x axis)
    deltaM = Maero + Mcontrol; % M moment sum (abt y axis)
    deltaN = Naero + Ncontrol; % N moment sum (abt Z axis)
    
    if t > 2 && t < 2.5
        deltaL = deltaL + Disturb(2); % Disturbed L moment (abt x axis)
        deltaM = deltaM + Disturb(1); % Disturbed M moment (abt y axis)
        deltaN = deltaN + Disturb(3); % Disturbed N moment (abt z axis)
    end
    Xaero = - eta^2 * deltau^2 * sign(deltau); % x component of aerodynamic force
    Yaero = - eta^2 * deltav^2 * sign(deltav); % y component of aerodynamic force
    Zaero = - eta^2 * deltaw^2 * sign(deltaw); % z component of aerodynamic force
    Xcontrol = 0; % x control
    Ycontrol = 0; % y control
    Zcontrol = -sum(Force); % gravitational force counteraction
    deltaX = Xaero + Xcontrol; % sum of x forces
    deltaY = Yaero + Ycontrol; % sum of y forces
    
    
    deltaZ = 0; % sum of z forces
    deltapdot = ((-k_1*deltap)/I_x) -((k_2*deltaphi)/I_x -((k_2*k_3)*(deltavr-deltav))/I_x  ); % roll rate derrivative
    deltaqdot = ((-k_4*deltaq)/I_y) -((k_5*deltatheta)/I_y -((k_5*k_6)*(deltaur-deltau))/I_y); % pitch rate derrivative
    
    
    deltardot = deltaN/I_z; % yaw rate derrivative
  
    dOmega_bdt = [deltapdot, deltaqdot, deltardot]';
    deltaudot = -g*deltatheta; % 4.7,1 acceleration in the x axis
    deltavdot = g*deltaphi; % acceleration in the y axis
    deltawdot = deltaZ/m; % acceleration in the z axis
    dVbdt = [deltaudot, deltavdot, deltawdot]';
    xdot =  deltau;
    ydot =  deltav;
    zdot = deltaw;
    dVEdt = [xdot, ydot, zdot]';
    phidot  = deltap; % bank roc (pg 104)
    thetadot = deltaq; % elevation roc (pg 104)
    psidot   = deltar; % azimuth roc (pg 104)
    dEuldt = [phidot, thetadot, psidot]';
    dydt = [dVbdt; dOmega_bdt; dEuldt; dVEdt];
    
    
    
end