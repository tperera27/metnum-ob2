function Runge_Kutta (intervals, tmax, dt, x0, v0)
  
    %Calculo aceleración
    a = @(u) (-((100*u) + (1000*(u**3)))); 
  
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
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
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
        resRK(1,1) = phiRK{1,1}*xRK(ts) + phiRK{1,2}*vRK(ts) + ((0.5)*((tao{1}*dt)**2)*(alphaRK{1,1}*a(xRK(ts))));
        %v_{t_s + tao1*dt}
        resRK(1,2) = phiRK{1,1}*vRK(ts) + phiRK{1,2}*a(xRK(ts));

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
    
    %plot(time, xRK, "r");
    %xlabel ("Tiempo", "fontsize", 20);
    %ylabel ("Posición", "fontsize", 20);
    %title ("Resorte elástico - Método de Runge-Kutta de tercer orden", "fontsize", 20);
    
    plot3(xRK, vRK, time);
    %plot(xRK, vRK);
    xlabel ("Posición", "fontsize", 20);
    ylabel ("Velocidad", "fontsize", 20);
    zlabel ("Tiempo", "fontsize", 20);
    title ("Resorte elástico - Método de Runge-Kutta de tercer orden", "fontsize", 20);
end


