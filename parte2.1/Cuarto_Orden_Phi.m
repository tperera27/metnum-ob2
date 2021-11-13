function Cuarto_Orden (intervals, tmax, dt, x0, v0)
    phi = cell(4,5);
    tao1 = 1/2;
    tao2 = 1;
    tao3 = 1;


    a = @(u) (-((100*u) + (1000*(u**3))));
    

    phi{1,1} = 1;
    phi{1,2} = dt / 3;

    phi{2,1} = 1;
    phi{2,2} = dt/8;
    phi{2,3} = (3/8)*dt;

    phi{3,1} = 1;
    phi{3,2} = (dt/2);
    phi{3,3} = -(3/2)*dt;
    phi{3,4} = 2*dt;
    
    phi{4,1} = 1;
    phi{4,2} = (dt/6);
    phi{4,3} = 0;
    phi{4,4} = (2/3)*dt;
    phi{4,5} = (dt/6);

    tao = {1/3, 1/2, 1};
    
    time = 0:dt:tmax;
    x = zeros(1, intervals);
    x(1) = x0;
    v = zeros(1, intervals);
    v(1) = v0;
    res4 = zeros(3,2);

    alpha = cell(4,4);
    alpha{1,1} = 1;
    alpha{2,1} = -3/5;
    alpha{2,2} = 3/5;
    alpha{3,1} = 3/5;
    alpha{3,2} = 3/5;
    alpha{3,3} = 0;
    alpha{4,1} = 0;
    alpha{4,2} = 0;
    alpha{4,3} = 0;
    alpha{4,4} = 0;
        
    ts = 1;
    while (ts < intervals)
        %Ecuacion 4 y 5
        %x_{t_s + dt/3}
        res4(1,1) = phi{1, 1}*x(ts) + phi{1,2}*v(ts) + (1/2)*((dt/3)^2)*((alpha{1,1})*a(x(ts)));
        %v_{t_s + dt/3}
        res4(1,2) = phi{1, 1}*v(ts) + phi{1,2}*a(x(ts));

        %Ecuacion 8 y 9
        %x_{t_s + dt/2}
        res4(2,1) = phi{2, 1}*x(ts) + phi{2,2}*v(ts) + phi{2,3}*res4(1,2) +  (1/2)*((dt/2)^2)*( (alpha{2,1}*a(x(ts))) + (alpha{2,2}*a(res4(1,1))) );
        %v_{t_s + dt/2}
        res4(2,2) = phi{2, 1}*v(ts) + phi{2,2}*a(x(ts)) + phi{2,3}*a(res4(1,1));
        
        %Ecuacion 28 y 29
        %x_{t_s + dt}
        res4(3,1) = phi{3, 1}*x(ts) + phi{3,2}*v(ts) + phi{3,3}*res4(1,2) + phi{3,4}*res4(2,2) + (1/2)*(dt^2)*( (alpha{3,1}*a(x(ts))) + (alpha{3,2}*a(res4(1,1))));
        %v_{t_s + dt}
        res4(3,2) = phi{3, 1}*v(ts) + phi{3,2}*a(x(ts)) + phi{3, 3}*a(res4(1,1)) + phi{3,4}*a(res4(2,1));
        
        %Ecuacion 30 y 31
        %x_{t_s + dt}
        x(ts + 1) = phi{4, 1}*x(ts) + phi{4,2}*v(ts) + phi{4,3}*res4(1,2) + phi{4,4}*res4(2,2) + phi{4,5}*res4(3,2);
        %v_{t_s + dt}
        v(ts + 1) = phi{4, 1}*v(ts) + phi{4,2}*a(x(ts)) + phi{4,3}*a(res4(1,1)) + phi{4,4}*a(res4(2,1)) + phi{4,5}*a(res4(3,1));
        
        ts = ts + 1;
    endwhile
    
    plot(time, x,"r");
    xlabel ("t");
    ylabel ("u(t)");
    title ("");
end


