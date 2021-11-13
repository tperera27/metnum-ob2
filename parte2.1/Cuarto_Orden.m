function Cuarto_Orden (intervals, tmax, dt, x0, v0)
    
    %Calculo aceleración
    a = @(u) (-((100*u) + (1000*(u**3))));
    
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
        resCO(1,1) = xCO(ts) + (dt/3)*vCO(ts) + (1/18)*a(xCO(ts));
        %v_{t_s + dt/3}
        resCO(1,2) = vCO(ts) + (dt/3)*a(xCO(ts));

        %Ecuacion 8 y 9
        %x_{t_s + dt/2}
        resCO(2,1) = xCO(ts) + (dt/2)*vCO(ts) + ((dt^2)/8)*( ((2/5)*a(xCO(ts))) + ((3/5)*a(resCO(1,1))) );
        %v_{t_s + dt/2}
        resCO(2,2) = vCO(ts) + (dt/8)*a(xCO(ts)) + ((3*dt)/8)*a(resCO(2,1));
        
        %Ecuacion 28 y 29
        %x_{t_s + dt}
        resCO(3,1) = xCO(ts) + dt*vCO(ts) + (1/2)*(dt^2)*( ((1/10)*a(xCO(ts))) - ((3/5)*a(resCO(1,1))) + ((3/2)*a(resCO(2,1))) );
        %v_{t_s + dt}
        resCO(3,2) = vCO(ts) + (1/2)*dt*a(xCO(ts)) + 2*dt*a(resCO(2,1));
        
        %Ecuacion 30 y 31
        %x_{t_s + dt}
        xCO(ts + 1) = xCO(ts) + dt*vCO(ts) + (1/2)*(dt^2)*( (2*a(xCO(ts))) - (3*a(resCO(1,1))) + (7*a(resCO(2,1))) );
        %v_{t_s + dt}
        vCO(ts + 1) = vCO(ts) + (1/6)*dt*a(xCO(ts)) + (2/3)*dt*a(resCO(2,1)) + (1/6)*dt*a(resCO(3,1));
        
        ts = ts + 1;
    endwhile
    
    plot(time, xCO,"r");
    xlabel ("t");
    ylabel ("u(t)");
    title ("");
end


