function Tercer_Orden (intervals, tmax, dt, x0, v0)
    
    %Calculo aceleración
    a = @(u) (-((100*u) + (1000*(u**3))));
    
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
    
    %Vector con los tiempos a graficar    
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
        resTO(1,1) = xTO(ts) + (dt/3)*vTO(ts) + ((0.5)*((dt^2)/9)*(a(xTO(ts))));
        %v_{t_s + tao1*dt}
        resTO(1,2) = vTO(ts) + (dt/3)*a(xTO(ts));

        %Segundo orden
        %x_{t_s + tao2*dt}
        resTO(2,1) = xTO(ts) + ((2*dt)/3)*vTO(ts) + ((dt^2)/27)*( (2*a(xTO(ts))) + (4*a(resTO(1,1))) );
        %v_{t_s + tao2*dt}
        resTO(2,2) = vTO(ts) + (2/3)*dt*a(resTO(1,1));
        
        %Tercer orden
        %x_{t_s + dt}
        xTO(ts + 1) = phiTO{3,1}*xTO(ts) + phiTO{3,2}*vTO(ts) + phiTO{3,3}*resTO(1,2) + phiTO{3,4}*resTO(2,2) + ( 0.5*(dt**2)*( (alphaTO{3,1}*a(xTO(ts))) + (alphaTO{3,2}*a(resTO(1,1))) +(alphaTO{3,3}*a(resTO(2,1))) ) );
        %v_{t_s + dt}
        vTO(ts + 1) = phiTO{3,1}*vTO(ts) + phiTO{3,2}*a(xTO(ts)) + phiTO{3,3}*a(resTO(1,1)) + phiTO{3,4}*a(resTO(2,1));
        
        ts = ts + 1;
    endwhile
    
    %plot(time, xTO, "r.-");
    %xlabel ("Tiempo", "fontsize", 20);
    %ylabel ("Posición", "fontsize", 20);
    %title ("Resorte elástico - Método explícito de tercer orden", "fontsize", 20);
    
    %plot3(xTO, vTO, time);
    plot(xTO, vTO);
    xlabel ("Posición", "fontsize", 20);
    ylabel ("Velocidad", "fontsize", 20);
    %zlabel ("Tiempo", "fontsize", 20);
    title ("Resorte elástico - Método explícito de tercer orden", "fontsize", 20);
end


