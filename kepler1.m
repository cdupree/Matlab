function kepler1


mu_earth = 11472;
mu_moon = 141.3;

t = 0:0.001:1;

% state is x y z v_x v_y v_z
% 42.5212 
init_state = [6.5917 0 0 0 42.5212 0];
init_state = [6.51 0 0 0 42 0];

r = [init_state(1) init_state(2) init_state(3)];
v = [init_state(4) init_state(5) init_state(6)];
h = cross(r,v);
e = cross(v,h)/mu_earth - r/norm(r);

options=odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',.001);
[t,state] = ode45(@rhs, t, init_state, options);

plot(state(:,1), state(:,2));
txt = {'eccentricty = ', norm(e)};
text(5, -6, txt);
hold on;
p = plot(state(1,1),state(1,2),'o','MarkerFaceColor','red');
axis manual;
tPos = text(-7,7, sprintf("t = %f", 0));
xPos = text(-7,6, sprintf("x = %f", state(1,1)));
yPos = text(-7,5, sprintf("y = %f", state(1,2)));
for ti = 2:length(state)
    p.XData = state(ti,1);
    p.YData = state(ti,2);
    set(tPos, 'String', sprintf("t = %f", t(ti)));
    set(xPos, 'String', sprintf("x = %f", state(ti,1)));
    set(yPos, 'String', sprintf("y = %f", state(ti,2)));
    drawnow;
    pause(.05);
end

    function state_dot = rhs(t, state)     
        r = sqrt(state(1)^2 + state(2)^2 + state(3)^2);
        f_mag = mu_earth / r^3;
        
        dxdt = state(4);
        dydt = state(5);
        dzdt = state(6);
        d2xdt2 =  -1 * f_mag * state(1);
        d2ydt2 =  -1 * f_mag * state(2);
        d2zdt2 =  -1 * f_mag * state(3);
        
        state_dot = [dxdt; dydt; dzdt; d2xdt2; d2ydt2; d2zdt2];
    end
end
