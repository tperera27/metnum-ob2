function Runge_Kutta (intervals, tmax, dt, x0, v0)
    phi = cell(3,4);
    tao1 = 1/2;
    tao2 = 1;


    a = @(u) (-((100*u) + (1000*(u^3))));
    

    phi{1,1} = @(t) (1);
    phi{1,2} = @(t) (t);

    phi{2,1} = @(t) (1);
    phi{2,2} = @(t) (-(t*(t + (tao1 * dt))) / (3 * tao1 * dt));
    phi{2,3} = @(t) ((t*t) / (tao1*dt));

    phi{3,1} = @(t) (1);
    phi{3,2} = @(t) ((t*(2*(t^2) - 3*tao1*dt*t -3*tao2*dt*t + 6*tao1*tao2*(dt^2))) / (6*tao1*tao2*(dt^2)));
    phi{3,3} = @(t) (-((t^2) - 2*t + 3*tao2*dt) / (6*tao1*(-tao2+tao1)*(dt^2)));
    phi{3,4} = @(t) ( ((t^2) - 2*t + 3*tao1*dt) / (6*tao2*(-tao2+tao1)*(dt^2)));

    tao = {1/2, 1, 1};
    
    time = 0:dt:tmax;
    x = zeros(1, intervals);
    x(1) = x0;
    v = zeros(1, intervals);
    v(1) = v0;
    resRK = zeros(2,2);

    alpha = cell(3,3);
    alpha{1,1} = 0;
    alpha{2,1} = 0;
    alpha{2,2} = 0;
    alpha{3,1} = 0;
    alpha{3,2} = 0;
    alpha{3,3} = 0;
  
    ts = 1;
    while (ts < intervals)
        %Primer orden
        %x_{t_s + tao1*dt}
        resRK(1,1) = phi{1, 1}(tao{1}*dt)*x(ts) + phi{1,2}(tao{1}*dt)*v(ts) + ((0.5)*(dt^2)*(alpha{1,1}*a(x(ts))));
        %v_{t_s + tao1*dt}
        resRK(1,2) = phi{1, 1}(tao{1}*dt)*v(ts) + phi{1,2}(tao{1}*dt)*a(x(ts));

        %Segundo orden
        %x_{t_s + tao2*dt}
        resRK(2,1) = phi{2,1}(dt*tao{2})*x(ts) + phi{2,2}(dt*tao{2})*v(ts) + phi{2,3}(dt*tao{2})*resRK(1,2) + ( 0.5*(dt^2)*( (alpha{2,1}*a(x(ts))) + (alpha{2,2}*a(resRK(1,1))) ) );
        %v_{t_s + tao2*dt}
        resRK(2,2) = phi{2,1}(tao{2}*dt)*v(ts) + phi{2,2}(tao{2}*dt)*a(x(ts)) + phi{2,3}(tao{2}*dt)*a(resRK(1,1));
        
        %Tercer orden
        %x_{t_s + dt}
        x(ts + 1) = phi{3,1}(dt)*x(ts) + phi{3,2}(dt)*v(ts) + phi{3,3}(dt)*resRK(1,2) + phi{3,4}(dt)*resRK(2,2) + ( 0.5*(dt^2)*( (alpha{3,1}*a(x(ts))) + (alpha{3,2}*a(resRK(1,1))) +(alpha{3,3}*a(resRK(2,1))) ) );
        %v_{t_s + dt}
        v(ts + 1) = phi{3,1}(dt)*x(ts) + phi{3,2}(dt)*v(ts) + phi{3,3}(dt)*resRK(1,2) + phi{3,4}(dt)*a(resRK(1,2));

        ts = ts + 1;
    endwhile
    
    display(time);
    display(x);
    plot(time, x,"r");
    xlabel ("ts");
    ylabel ("u(ts)");
    title ("Simple 2-D Plot");
end


