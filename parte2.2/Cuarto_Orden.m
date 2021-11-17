function Cuarto_Orden (intervals, tmax, dt, x0, v0)   
    
    %Calculo aceleración
    a = @(tita) (-sin(tita));
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Posición angular
    xCO = zeros(1, intervals);
    xCO(1) = x0;
    
    %Velocidad angular
    vCO = zeros(1, intervals);
    vCO(1) = v0;
    
    %Resultados intermedios
    resCO = zeros(3,2);
        
    ts = 1;
    while (ts < intervals)
        %Ecuacion 4 y 5
        %x_{t_s + dt/3}
        resCO(1,1) = xCO(ts) + (dt/3)*vCO(ts) + (1/18)*(dt^2)*a(xCO(ts));
        %v_{t_s + dt/3}
        resCO(1,2) = vCO(ts) + (dt/3)*a(xCO(ts));

        %Ecuacion 8 y 9
        %x_{t_s + dt/2}
        resCO(2,1) = xCO(ts) + (dt/8)*vCO(ts) + ((3/8)*dt)*resCO(1,2) + ((dt^2)/8)*(- ((3/5)*a(xCO(ts))) + ((3/5)*a(resCO(1,1))) );
        %v_{t_s + dt/2}
        resCO(2,2) = vCO(ts) + (dt/8)*a(xCO(ts)) + ((3*dt)/8)*a(resCO(1,1));
        
        %Ecuacion 28 y 29
        %x_{t_s + dt}
        resCO(3,1) = xCO(ts) + (1/2)*dt*vCO(ts) - (3/2)*dt*resCO(1,2) + 2*dt*resCO(2,2) + (1/10)*(dt^2)*( (3*a(xCO(ts))) - (3*a(resCO(1,1))) );
        %v_{t_s + dt}
        resCO(3,2) = vCO(ts) + (1/2)*dt*a(xCO(ts)) - (3/2)*dt*a(resCO(1,1)) + 2*dt*a(resCO(2,1));
        
        %Ecuacion 30 y 31
        %x_{t_s + dt}
        xCO(ts + 1) = xCO(ts) + (dt/6)*vCO(ts) + (2/3)*dt*resCO(2,2) + (1/6)*dt*resCO(3,2);
        %v_{t_s + dt}
        vCO(ts + 1) = vCO(ts) + (1/6)*dt*a(xCO(ts)) + (2/3)*dt*a(resCO(2,1)) + (1/6)*dt*a(resCO(3,1));
        
        ts = ts + 1;
    endwhile
    
    plot(time, xCO,"r");
    xlabel ("t");
    ylabel ("u(t)");
    title ("");
end


