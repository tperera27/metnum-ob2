function Graficas (intervals, tmax, dt, x0, v0)
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Calculo aceleración
    a = @(u) (-((100*u) + (1000*(u**3)))); 
    
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
    %Phi precalculado
    phiCO = cell(4,5);    
    
    phiCO{1,1} = 1;
    phiCO{1,2} = dt / 3;

    phiCO{2,1} = 1;
    phiCO{2,2} = dt/8;
    phiCO{2,3} = (3/8)*dt;

    phiCO{3,1} = 1;
    phiCO{3,2} = (dt/2);
    phiCO{3,3} = -(3/2)*dt;
    phiCO{3,4} = 2*dt;
    
    phiCO{4,1} = 1;
    phiCO{4,2} = (dt/6);
    phiCO{4,3} = 0;
    phiCO{4,4} = (2/3)*dt;
    phiCO{4,5} = (dt/6);

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
        resCO(1,1) = phiCO{1, 1}*xCO(ts) + phiCO{1,2}*vCO(ts) + (1/2)*((dt/3)^2)*((alphaCO{1,1})*a(xCO(ts)));
        %v_{t_s + dt/3}
        resCO(1,2) = phiCO{1, 1}*vCO(ts) + phiCO{1,2}*a(xCO(ts));

        %Ecuacion 8 y 9
        %x_{t_s + dt/2}
        resCO(2,1) = phiCO{2, 1}*xCO(ts) + phiCO{2,2}*vCO(ts) + phiCO{2,3}*resCO(1,2) +  (1/2)*((dt/2)^2)*( (alphaCO{2,1}*a(xCO(ts))) + (alphaCO{2,2}*a(resCO(1,1))) );
        %v_{t_s + dt/2}
        resCO(2,2) = phiCO{2, 1}*vCO(ts) + phiCO{2,2}*a(xCO(ts)) + phiCO{2,3}*a(resCO(1,1));
        
        %Ecuacion 28 y 29
        %x_{t_s + dt}
        resCO(3,1) = phiCO{3, 1}*xCO(ts) + phiCO{3,2}*vCO(ts) + phiCO{3,3}*resCO(1,2) + phiCO{3,4}*resCO(2,2) + (1/2)*(dt^2)*( (alphaCO{3,1}*a(xCO(ts))) + (alphaCO{3,2}*a(resCO(1,1))));
        %v_{t_s + dt}
        resCO(3,2) = phiCO{3, 1}*vCO(ts) + phiCO{3,2}*a(xCO(ts)) + phiCO{3, 3}*a(resCO(1,1)) + phiCO{3,4}*a(resCO(2,1));
        
        %Ecuacion 30 y 31
        %x_{t_s + dt}
        xCO(ts + 1) = phiCO{4, 1}*xCO(ts) + phiCO{4,2}*vCO(ts) + phiCO{4,3}*resCO(1,2) + phiCO{4,4}*resCO(2,2) + phiCO{4,5}*resCO(3,2);
        %v_{t_s + dt}
        vCO(ts + 1) = phiCO{4, 1}*vCO(ts) + phiCO{4,2}*a(xCO(ts)) + phiCO{4,3}*a(resCO(1,1)) + phiCO{4,4}*a(resCO(2,1)) + phiCO{4,5}*a(resCO(3,1));
        
        ts = ts + 1;
    endwhile
    
    plot(time, xRK, "r", time, xTO, "b", time, xCO, "g")
    xlabel ("Tiempo t", "fontsize", 20);
    ylabel ("Posición", "fontsize", 20);
    h = legend ("Runge-Kutta", "Tercer orden", "Cuarto orden");
    legend (h, "location", "northeastoutside");
    set (h, "fontsize", 20)
    title ("Resorte elástico", "fontsize", 30);

 
end


