function Cuarto_Orden_Phi (intervals, tmax, dt, x0, v0)
    
    %Calculo aceleración
    a = @(tita) (-sin(tita));
    
    %Phi precalculado
    phi = cell(4,5);
    
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
    
    %Parámetros  
    tao = {1/3, 1/2, 1};
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
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Posición angular
    x = zeros(1, intervals);
    x(1) = x0;
    
    %Velocidad angular
    v = zeros(1, intervals);
    v(1) = v0;
    
    %Resultados intermedios
    resCO = zeros(3,2);
    
        
    ts = 1;
    while (ts < intervals)
        %Ecuacion 4 y 5
        %x_{t_s + dt/3}
        resCO(1,1) = phi{1, 1}*x(ts) + phi{1,2}*v(ts) + (1/2)*((dt/3)^2)*((alpha{1,1})*a(x(ts)));
        %v_{t_s + dt/3}
        resCO(1,2) = phi{1, 1}*v(ts) + phi{1,2}*a(x(ts));

        %Ecuacion 8 y 9
        %x_{t_s + dt/2}
        resCO(2,1) = phi{2, 1}*x(ts) + phi{2,2}*v(ts) + phi{2,3}*resCO(1,2) +  (1/2)*((dt/2)^2)*( (alpha{2,1}*a(x(ts))) + (alpha{2,2}*a(resCO(1,1))) );
        %v_{t_s + dt/2}
        resCO(2,2) = phi{2, 1}*v(ts) + phi{2,2}*a(x(ts)) + phi{2,3}*a(resCO(1,1));
        
        %Ecuacion 28 y 29
        %x_{t_s + dt}
        resCO(3,1) = phi{3, 1}*x(ts) + phi{3,2}*v(ts) + phi{3,3}*resCO(1,2) + phi{3,4}*resCO(2,2) + (1/2)*(dt^2)*( (alpha{3,1}*a(x(ts))) + (alpha{3,2}*a(resCO(1,1))));
        %v_{t_s + dt}
        resCO(3,2) = phi{3, 1}*v(ts) + phi{3,2}*a(x(ts)) + phi{3, 3}*a(resCO(1,1)) + phi{3,4}*a(resCO(2,1));
        
        %Ecuacion 30 y 31
        %x_{t_s + dt}
        x(ts + 1) = phi{4, 1}*x(ts) + phi{4,2}*v(ts) + phi{4,3}*resCO(1,2) + phi{4,4}*resCO(2,2) + phi{4,5}*resCO(3,2);
        %v_{t_s + dt}
        v(ts + 1) = phi{4, 1}*v(ts) + phi{4,2}*a(x(ts)) + phi{4,3}*a(resCO(1,1)) + phi{4,4}*a(resCO(2,1)) + phi{4,5}*a(resCO(3,1));
        
        ts = ts + 1;
    endwhile
    
    plot(time, x, "g");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo simple - Método explicito de cuarto orden", "fontsize", 30);
end


