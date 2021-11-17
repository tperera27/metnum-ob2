function Graficas (intervals, tmax, dt, x0, v0)
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Calculo aceleración
    a = @(tita) (-sin(tita)); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNGE KUTTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %Phi precalculado
    phiRK = cell(3,4);
    
    phiRK{1,1} = 1;
    phiRK{1,2} = dt / 2;

    phiRK{2,1} = 1;
    phiRK{2,2} = -(dt);
    phiRK{2,3} = 2*dt;

    phiRK{3,1} = 1;
    phiRK{3,2} = (dt/6);
    phiRK{3,3} = (2/3)*dt;
    phiRK{3,4} = (dt/6);

    %Parámetros
    alphaRK = cell(3,3);
    alphaRK{1,1} = 0;
    alphaRK{2,1} = 0;
    alphaRK{2,2} = 0;
    alphaRK{3,1} = 0;
    alphaRK{3,2} = 0;
    alphaRK{3,3} = 0;
    tao = {1/2, 1, 1};
    
    %Posición
    xRK = zeros(1, intervals);
    xRK(1) = x0;
    
    %Velocidad
    vRK = zeros(1, intervals);
    vRK(1) = v0;
    
    %Resultados intermedios
    resRK = zeros(2,2);
        
    ts = 1;
    while (ts < intervals)
        %Primer orden
        %x_{t_s + tao1*dt}
        resRK(1,1) = phiRK{1, 1}*xRK(ts) + phiRK{1,2}*vRK(ts) + ((0.5)*((tao{1}*dt)**2)*(alphaRK{1,1}*a(xRK(ts))));
        %v_{t_s + tao1*dt}
        resRK(1,2) = phiRK{1, 1}*vRK(ts) + phiRK{1,2}*a(xRK(ts));

        %Segundo orden
        %x_{t_s + tao2*dt}
        resRK(2,1) = phiRK{2,1}*xRK(ts) + phiRK{2,2}*vRK(ts) + phiRK{2,3}*resRK(1,2) + ( 0.5*((tao{2}*dt)**2)*( (alphaRK{2,1}*a(xRK(ts))) + (alphaRK{2,2}*a(resRK(1,1))) ) );
        %v_{t_s + tao2*dt}
        resRK(2,2) = phiRK{2,1}*vRK(ts) + phiRK{2,2}*a(xRK(ts)) + phiRK{2,3}*a(resRK(1,1));
        
        %Tercer orden
        %x_{t_s + dt}
        xRK(ts + 1) = phiRK{3,1}*xRK(ts) + phiRK{3,2}*vRK(ts) + phiRK{3,3}*resRK(1,2) + phiRK{3,4}*resRK(2,2) + ( 0.5*(dt**2)*( (alphaRK{3,1}*a(xRK(ts))) + (alphaRK{3,2}*a(resRK(1,1))) +(alphaRK{3,3}*a(resRK(2,1))) ) );
        %v_{t_s + dt}
        vRK(ts + 1) = phiRK{3,1}*vRK(ts) + phiRK{3,2}*a(xRK(ts)) + phiRK{3,3}*a(resRK(1,1)) + phiRK{3,4}*a(resRK(2,1));
        
        ts = ts + 1;
    endwhile
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TERCER ORDEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Phi precalculado
    phiTO = cell(3,4);
    
    phiTO{1,1} = 1;
    phiTO{1,2} = dt/3;

    phiTO{2,1} = 1;
    phiTO{2,2} = 0;
    phiTO{2,3} = (2*dt) / 3;

    phiTO{3,1} = 1;
    phiTO{3,2} = (dt/4);
    phiTO{3,3} = 0;
    phiTO{3,4} = ((3*dt)/4);

    %Parámetros    
    taoTO = {1/3, 2/3, 1};
    alphaTO = cell(3,3);
    alphaTO{1,1} = 1;
    alphaTO{2,1} = -2/3;
    alphaTO{2,2} = 2/3;
    alphaTO{3,1} = 1/3;
    alphaTO{3,2} = -2/3;
    alphaTO{3,3} = 1/3;
    
    %Vecor con los tiempos a graficar    
    time = 0:dt:tmax;
    
    %Posición    
    xTO = zeros(1, intervals);
    xTO(1) = x0;
    
    %Velocidad
    vTO = zeros(1, intervals);
    vTO(1) = v0;
    
    %Resultados intermedios
    resTO = zeros(2,2);
  
    ts = 1;
    while (ts < intervals)
        %Primer orden
        %x_{t_s + tao1*dt}
        resTO(1,1) = phiTO{1, 1}*xTO(ts) + phiTO{1,2}*vTO(ts) + ((0.5)*((taoTO{1}*dt)**2)*(alphaTO{1,1}*a(xTO(ts))));
        %v_{t_s + tao1*dt}
        resTO(1,2) = phiTO{1, 1}*vTO(ts) + phiTO{1,2}*a(xTO(ts));

        %Segundo orden
        %x_{t_s + tao2*dt}
        resTO(2,1) = phiTO{2,1}*xTO(ts) + phiTO{2,2}*vTO(ts) + phiTO{2,3}*resTO(1,2) + ( 0.5*((taoTO{2}*dt)**2)*( (alphaTO{2,1}*a(xTO(ts))) + (alphaTO{2,2}*a(resTO(1,1))) ) );
        %v_{t_s + tao2*dt}
        resTO(2,2) = phiTO{2,1}*vTO(ts) + phiTO{2,2}*a(xTO(ts)) + phiTO{2,3}*a(resTO(1,1));
        
        %Tercer orden
        %x_{t_s + dt}
        xTO(ts + 1) = phiTO{3,1}*xTO(ts) + phiTO{3,2}*vTO(ts) + phiTO{3,3}*resTO(1,2) + phiTO{3,4}*resTO(2,2) + ( 0.5*(dt**2)*( (alphaTO{3,1}*a(xTO(ts))) + (alphaTO{3,2}*a(resTO(1,1))) +(alphaTO{3,3}*a(resTO(2,1))) ) );
        %v_{t_s + dt}
        vTO(ts + 1) = phiTO{3,1}*vTO(ts) + phiTO{3,2}*a(xTO(ts)) + phiTO{3,3}*a(resTO(1,1)) + phiTO{3,4}*a(resTO(2,1));
        
        ts = ts + 1;
    endwhile
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUARTO ORDEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Posición
    xCO = zeros(1, intervals);
    xCO(1) = x0;
    
    %Velocidad
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
    
    plot(time, xRK, "k.-", time, xTO, "ro-", time, xCO, "b.-") 
    xlabel ("Tiempo t", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    h = legend ("Runge-Kutta", "Tercer orden", "Cuarto orden");
    legend (h, "location", "northeastoutside");
    set (h, "fontsize", 20)
    title ("Péndulo simple", "fontsize", 30);

 
end


