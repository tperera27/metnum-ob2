function Cuarto_Orden (intervals, tmax, dt, x0, v0)
    
    %Calculo aceleración
    a = @(u) (-((100*u) + (1000*(u**3))));
    

    %Parámetros    
    taoCO = {1/3, 1/2, 1};
    alphaCO = cell(4,4);
    alphaCO{1,1} = 1;
    alphaCO{2,1} = -3/5;
    alphaCO{2,2} = 3/5;
    alphaCO{3,1} = 3/5;
    alphaCO{3,2} = 3/5;
    alphaCO{3,3} = 0;
    alphaCO{4,1} = 0;
    alphaCO{4,2} = 0;
    alphaCO{4,3} = 0;
    alphaCO{4,4} = 0;
    
    %Vecor con los tiempos a graficar     
    time = 0:dt:tmax;
    
    %Posición
    xCO = zeros(1, intervals);
    xCO(1) = x0;
    
    %Velocidad
    vCO = zeros(1, intervals);
    vCO(1) = v0;
    
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
    
    
    %plot(time, xCO, "g");
    %xlabel ("Tiempo", "fontsize", 20);
    %ylabel ("Posición", "fontsize", 20);
    %title ("Resorte elástico - Método explícito de cuarto orden", "fontsize", 30);
    
    plot3(xCO, vCO, time);
    %plot(xCO, vCO);
    xlabel ("Posición", "fontsize", 20);
    ylabel ("Velocidad", "fontsize", 20);
    zlabel ("Tiempo", "fontsize", 20);
    title ("Resorte elástico - Método explícito de cuarto orden", "fontsize", 30);
end


