function Tercer_Orden(intervals, tmax, dt, x0, v0)
    
    %Calculo aceleración
    a = @(tita) (-sin(tita));
    
    %Phi precalculado
    phi = cell(3,4);
    
    phi{1,1} = 1;
    phi{1,2} = dt / 3;

    phi{2,1} = 1;
    phi{2,2} = 0;
    phi{2,3} = (2*dt) / 3;

    phi{3,1} = 1;
    phi{3,2} = (dt/4);
    phi{3,3} = 0;
    phi{3,4} = ((3*dt)/4);
    
    %Parámetros
    tao = {1/3, 2/3, 1};
    alpha = cell(3,3);
    alpha{1,1} = 1;
    alpha{2,1} = -2/3;
    alpha{2,2} = 2/3;
    alpha{3,1} = 1/3;
    alpha{3,2} = -2/3;
    alpha{3,3} = 1/3;
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Posición angular
    x = zeros(1, intervals);
    x(1) = x0;
    
    %Velocidad angular
    v = zeros(1, intervals);
    v(1) = v0;
    
    %Resultados intermedios
    resTO3 = zeros(2,2);

 
    ts = 1;
    while (ts < intervals)
        %Primer orden
        %x_{t_s + dt/3}
        resTO3(1,1) = phi{1, 1}*x(ts) + phi{1,2}*v(ts) + ((0.5)*((tao{1}*dt)**2)*(alpha{1,1}*a(x(ts))));
        %v_{t_s + tao1*dt}
        resTO3(1,2) = phi{1, 1}*v(ts) + phi{1,2}*a(x(ts));

        %Segundo orden
        %x_{t_s + 2dt/3}
        resTO3(2,1) = phi{2,1}*x(ts) + phi{2,2}*v(ts) + phi{2,3}*resTO3(1,2) + ( 0.5*((tao{2}*dt)**2)*( (alpha{2,1}*a(x(ts))) + (alpha{2,2}*a(resTO3(1,1))) ) );
        %v_{t_s + tao2*dt}
        resTO3(2,2) = phi{2,1}*v(ts) + phi{2,2}*a(x(ts)) + phi{2,3}*a(resTO3(1,1));
        
        %Tercer orden
        %x_{t_s + dt}
        x(ts + 1) = phi{3,1}*x(ts) + phi{3,2}*v(ts) + phi{3,3}*resTO3(1,2) + phi{3,4}*resTO3(2,2) + ( 0.5*(dt**2)*( (alpha{3,1}*a(x(ts))) + (alpha{3,2}*a(resTO3(1,1))) +(alpha{3,3}*a(resTO3(2,1))) ) );
        %v_{t_s + dt}
        v(ts + 1) = phi{3,1}*v(ts) + phi{3,2}*a(x(ts), v(ts)) + phi{3,3}*a(resTO3(1,1)) + phi{3,4}*a(resTO3(2,1));
        
        ts = ts + 1;
    endwhile
    
    plot(time, x, "b");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo simple - Método explicito de tercer orden", "fontsize", 30);
end


