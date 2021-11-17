function Graficas (intervals, tmax, dt, r0, v0, tita0, w0)
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Calculo aceleración del resorte (radial)
    art = @(rt,vt,titat,wt) (((wt^2)*(0.5 + rt) + (cos(titat))*9.81) - (98.1 * rt));
    
    a = @(tita, w, r) ((w^2)*(0.5 + r) + cos(tita)*9.81 - 98.1 * r);
    
    %Calculo aceleración del pendulo (angular)
    aat = @(rt,vt,titat,wt) (-( (2*wt*vt) + (9.81*sin(titat)) ) / (0.5 + rt) );
    
    tita_a = @(tita, w, r, v) ((-2*v*w - 9.81*sin(tita)) / (0.5 + r));
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNGE KUTTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %phiRK precalculado
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
    
    %Posición radial
    rRK = zeros(1, intervals);
    rRK(1) = r0;
    
    %Velocidad radial
    vRK = zeros(1, intervals);
    vRK(1) = v0;
    
    %Posición angular
    titaRK = zeros(1, intervals);
    titaRK(1) = tita0;
    
    %Velocidad angular
    wRK = zeros(1, intervals);
    wRK (1) = w0;
    
    %Resultados intermedios
    resResorteRK  = zeros(2,2);
    resPenduloRK = zeros(2,2);

        
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        % Resorte primer paso
        resResorteRK (1,1) = phiRK{1, 1}*rRK(ts) + phiRK{1,2}*vRK(ts);
        resResorteRK (1,2) = phiRK{1, 1}*vRK(ts) + phiRK{1,2}*a(titaRK(ts), wRK (ts), rRK(ts));

        % Pendulo primer paso
        resPenduloRK(1,1) = phiRK{1, 1}*titaRK(ts) + phiRK{1,2}*wRK (ts);
        resPenduloRK(1,2) = phiRK{1, 1}*wRK (ts) + phiRK{1,2}*tita_a(titaRK(ts), wRK (ts), rRK(ts), vRK(ts));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        %Resorte segundo paso
        resResorteRK (2,1) = phiRK{2,1}*rRK(ts) + phiRK{2,2}*vRK(ts) + phiRK{2,3}*resResorteRK (1,2);
        resResorteRK (2,2) = phiRK{2,1}*vRK(ts) + phiRK{2,2}*a(titaRK(ts), wRK (ts), rRK(ts)) + phiRK{2,3}*a(resPenduloRK(1,1), resPenduloRK(1,2), resResorteRK (1,1));
        
        %Pendulo segundo paso
        resPenduloRK(2,1) = phiRK{2,1}*titaRK(ts) + phiRK{2,2}*wRK (ts) + phiRK{2,3}*resPenduloRK(1,2);
        resPenduloRK(2,2) = phiRK{2,1}*wRK (ts) + phiRK{2,2}*tita_a(titaRK(ts), wRK (ts), rRK(ts), vRK(ts)) + phiRK{2,3}*tita_a(resPenduloRK(1,1), resPenduloRK(1,2), resResorteRK (1,1), resResorteRK (1,2));
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Resorte tercer paso
        rRK(ts + 1) = phiRK{3,1}*rRK(ts) + phiRK{3,2}*vRK(ts) + phiRK{3,3}*resResorteRK (1,2) + phiRK{3,4}*resResorteRK (2,2);
        vRK(ts + 1) = phiRK{3,1}*vRK(ts) + phiRK{3,2}*a(titaRK(ts), wRK (ts), rRK(ts)) + phiRK{3,3}*a(resPenduloRK(1,1), resPenduloRK(1,2), resResorteRK (1,1)) + phiRK{3,4}*a(resPenduloRK(2,1), resPenduloRK(2,2), resResorteRK (2,1));

        % Pendulo tercer paso
        titaRK(ts + 1) = phiRK{3,1}*titaRK(ts) + phiRK{3,2}*wRK (ts) + phiRK{3,3}*resPenduloRK(1,2) + phiRK{3,4}*resPenduloRK(2,2);
        wRK (ts + 1) = phiRK{3,1}*wRK (ts) + phiRK{3,2}*tita_a(titaRK(ts), wRK (ts), rRK(ts), vRK(ts)) + phiRK{3,3}*tita_a(resPenduloRK(1,1), resPenduloRK(1,2), resResorteRK (1,1), resResorteRK (1,2)) + phiRK{3,4}*tita_a(resPenduloRK(2,1), resPenduloRK(2,2), resResorteRK (2,1), resResorteRK (2,2));

        ts = ts + 1;
    endwhile
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TERCER ORDEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %phiTO precalculado
    phiTO = cell(3,4);
    
    phiTO{1,1} = 1;
    phiTO{1,2} = dt / 3;

    phiTO{2,1} = 1;
    phiTO{2,2} = 0;
    phiTO{2,3} = (2*dt) / 3;

    phiTO{3,1} = 1;
    phiTO{3,2} = (dt/4);
    phiTO{3,3} = 0;
    phiTO{3,4} = ((3*dt)/4);
    
    %Parámetros
    taoTO = {1/3, 2/3, 1};
    alphaTO = zeros(3,3);
    alphaTO(1,1) = 1;
    alphaTO(2,1) = -2/3;
    alphaTO(2,2) = 2/3;
    alphaTO(3,1) = 1/3;
    alphaTO(3,2) = -2/3;
    alphaTO(3,3) = 1/3;
    
    %Posición radial
    rTO = zeros(1, intervals);
    rTO(1) = r0;
    
    %Velocidad radial
    vTO = zeros(1, intervals);
    vTO(1) = v0;
    
    %Posición angular
    titaTO = zeros(1, intervals);
    titaTO(1) = tita0;
    
    %Velocidad angular
    wTO = zeros(1, intervals);
    wTO(1) = w0;
    
    %Resultados intermedios
    resResorteTO = zeros(2,2);
    resPenduloTO = zeros(2,2);
    
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 4 y 5
        
        % Resorte primer paso
        % u{ts + dt/3}
        resResorteTO(1,1) = phiTO{1,1}*rTO(ts) + phiTO{1,2}*vTO(ts) + (1/2)*((dt/3)^2)*(alphaTO(1,1))*art(rTO(ts),vTO(ts),titaTO(ts),wTO(ts));
        % v{ts + dt/3}
        resResorteTO(1,2) = phiTO{1,1}*vTO(ts) + phiTO{1,2}*art(rTO(ts),vTO(ts),titaTO(ts),wTO(ts));
        
        % Pendulo primer paso
        % u{ts + dt/3}
        resPenduloTO(1,1) = phiTO{1,1}*titaTO(ts) + phiTO{1,2}*wTO(ts) + (1/2)*((dt/3)^2)*(alphaTO(1,1))*aat(rTO(ts),vTO(ts),titaTO(ts),wTO(ts));
        % v{ts + dt/3}
        resPenduloTO(1,2) = phiTO{1,1}*wTO(ts) + phiTO{1,2}*aat(rTO(ts),vTO(ts),titaTO(ts),wTO(ts));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 8 y 9
        
        % Resorte segundo paso
        % u{ts + 2dt/3}
        resResorteTO(2,1) = phiTO{2,1}*rTO(ts) + phiTO{2,2}*vTO(ts) + phiTO{2,3}*resResorteTO(1,2) + (0.5)*((2*dt/3)^2)*( (alphaTO(2,1)*art(rTO(ts),vTO(ts),titaTO(ts),wTO(ts))) + (alphaTO(2,2)*art(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2))) );
        % v{ts + 2dt/3}
        resResorteTO(2,2) = phiTO{2,1}*vTO(ts) + phiTO{2,2}*art(rTO(ts),vTO(ts),titaTO(ts),wTO(ts)) + phiTO{2,3}*art(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2));

        % Pendulo segundo paso
        % u{ts + 2dt/3}
        resPenduloTO(2,1) = phiTO{2,1}*titaTO(ts) + phiTO{2,2}*wTO(ts) + phiTO{2,3}*resPenduloTO(1,2) + (0.5)*((2*dt/3)^2)*( (alphaTO(2,1)*aat(rTO(ts),vTO(ts),titaTO(ts),wTO(ts))) + (alphaTO(2,2)*aat(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2))) );
        % v{ts + 2dt/3}
        resPenduloTO(2,2) = phiTO{2,1}*wTO(ts) + phiTO{2,2}*aat(rTO(ts),vTO(ts),titaTO(ts),wTO(ts)) + phiTO{2,3}*aat(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 13 y 14
        
        % Resorte segundo paso
        % u{ts + dt}
        rTO(ts + 1) = phiTO{3,1}*rTO(ts) + phiTO{3,2}*vTO(ts) + phiTO{3,3}*resResorteTO(1,2) + phiTO{3,4}*resResorteTO(2,2) + (0.5)*(dt^2)*( (alphaTO(3,1)*art(rTO(ts),vTO(ts),titaTO(ts),wTO(ts))) + (alphaTO(3,2)*art(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2))) + (alphaTO(3,3)*art(resResorteTO(2,1), resResorteTO(2,2), resPenduloTO(2,1), resPenduloTO(2,2))) );
        % v{ts + dt}
        vTO(ts + 1) = phiTO{3,1}*vTO(ts) + phiTO{3,2}*art(rTO(ts),vTO(ts),titaTO(ts),wTO(ts)) + phiTO{3,3}*art(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2)) + phiTO{3,4}*art(resResorteTO(2,1), resResorteTO(2,2), resPenduloTO(2,1), resPenduloTO(2,2));
        
        % Pendulo segundo paso
        % u{ts + dt}
        titaTO(ts + 1) = phiTO{3,1}*titaTO(ts) + phiTO{3,2}*wTO(ts) + phiTO{3,3}*resPenduloTO(1,2) + phiTO{3,4}*resPenduloTO(2,2) + (0.5)*(dt^2)*( (alphaTO(3,1)*aat(rTO(ts),vTO(ts),titaTO(ts),wTO(ts))) + (alphaTO(3,2)*aat(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2))) + (alphaTO(3,3)*aat(resResorteTO(2,1), resResorteTO(2,2), resPenduloTO(2,1), resPenduloTO(2,2))) );
        % v{ts + dt}
        wTO(ts + 1) = phiTO{3,1}*wTO(ts) + phiTO{3,2}*aat(rTO(ts),vTO(ts),titaTO(ts),wTO(ts)) + phiTO{3,3}*aat(resResorteTO(1,1), resResorteTO(1,2), resPenduloTO(1,1), resPenduloTO(1,2)) + phiTO{3,4}*aat(resResorteTO(2,1), resResorteTO(2,2), resPenduloTO(2,1), resPenduloTO(2,2));
        
        ts = ts + 1;
    endwhile
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUARTO ORDEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Posición radial
    rCO = zeros(1, intervals);
    rCO(1) = r0;
    
    %Velocidad radial
    vCO = zeros(1, intervals);
    vCO(1) = v0;
    
    %Posición angular
    titaCO = zeros(1, intervals);
    titaCO(1) = tita0;
    
    %Velocidad angular
    wCO = zeros(1, intervals);
    wCO(1) = w0;
    
    %Resultados intermedios
    resResorteCO = zeros(3,2);
    resPenduloCO = zeros(3,2);
    
    %Calculo aceleración del resorte (radial)
    art = @(rt,vt,titat,wt) (((wt^2)*(0.5 + rt) + (cos(titat))*9.81) - (98.1 * rt));
    
    %Calculo aceleración del pendulo (angular)
    aat = @(rt,vt,titat,wt) (-( (2*wt*vt) + (9.81*sin(titat)) ) / (0.5 + rt) );
   
   
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 4 y 5
        
        % Resorte primer paso
        % u{ts + dt/3}
        resResorteCO(1,1) = rCO(ts) + (dt/3)*vCO(ts) + (1/18)*(dt^2)*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts));
        % v{ts + dt/3}
        resResorteCO(1,2) = vCO(ts) + (dt/3)*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts));
        
        % Pendulo primer paso
        % u{ts + dt/3}
        resPenduloCO(1,1) = titaCO(ts) + (dt/3)*wCO(ts) + (1/18)*(dt^2)*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts));
        % v{ts + dt/3}
        resPenduloCO(1,2) = wCO(ts) + (dt/3)*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 8 y 9
        
        % Resorte segundo paso
        % u{ts + dt/2}
        resResorteCO(2,1) = rCO(ts) + (dt/8)*vCO(ts) + ((3/8)*dt)*resResorteCO(1,2) + ((dt^2)/8)*(- ((3/5)*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts))) + ((3/5)*art(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2))) );
        % v{ts + dt/2}
        resResorteCO(2,2) = vCO(ts) + (dt/8)*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts)) + ((3*dt)/8)*art(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2));

        % Pendulo segundo paso
        % u{ts + dt/2}
        resPenduloCO(2,1) = titaCO(ts) + (dt/8)*wCO(ts) + ((3/8)*dt)*resPenduloCO(1,2) + ((dt^2)/8)*(- ((3/5)*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts))) + ((3/5)*aat(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2)))  );
        % v{ts + dt/2}
        resPenduloCO(2,2) = wCO(ts) + (dt/8)*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts)) + ((3*dt)/8)*aat(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuacion 28 y 29
        
        % Resorte tercer paso
        % u{ts + dt}
        resResorteCO(3,1) = rCO(ts) + (1/2)*dt*vCO(ts) - (3/2)*dt*resResorteCO(1,2) + 2*dt*resResorteCO(2,2) + (1/10)*(dt^2)*( (3*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts))) - (3*art(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2))) );
        % v{ts + dt}
        resResorteCO(3,2) = vCO(ts) + (1/2)*dt*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts)) - (3/2)*dt*art(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2)) + 2*dt*art(resResorteCO(2,1), resResorteCO(2,2), resPenduloCO(2,1), resPenduloCO(2,2));

        % Pendulo tercer paso
        % u{ts + dt}
        resPenduloCO(3,1) = titaCO(ts) + (1/2)*dt*wCO(ts) - (3/2)*dt*resPenduloCO(1,2) + 2*dt*resPenduloCO(2,2) + (1/10)*(dt^2)*( (3*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts))) - (3*aat(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2))) );
        % v{ts + dt}
        resPenduloCO(3,2) = wCO(ts) + (1/2)*dt*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts)) - (3/2)*dt*aat(resResorteCO(1,1), resResorteCO(1,2), resPenduloCO(1,1), resPenduloCO(1,2)) + 2*dt*aat(resResorteCO(2,1), resResorteCO(2,2), resPenduloCO(2,1), resPenduloCO(2,2));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Cuarto paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuacion 30 y 31
        
        % Resorte cuarto paso
        % u{ts + dt}
        rCO(ts + 1) = rCO(ts) + (dt/6)*vCO(ts) + (2/3)*dt*resResorteCO(2,2) + (1/6)*dt*resResorteCO(3,2);
        % v{ts + dt}
        vCO(ts + 1) = vCO(ts) + (1/6)*dt*art(rCO(ts),vCO(ts),titaCO(ts),wCO(ts)) + (2/3)*dt*art(resResorteCO(2,1), resResorteCO(2,2), resPenduloCO(2,1), resPenduloCO(2,2)) + (1/6)*dt*art(resResorteCO(3,1), resResorteCO(3,2), resPenduloCO(3,1), resPenduloCO(3,2));

        % Pendulo cuarto paso
        % u{ts + dt}
        titaCO(ts + 1) = titaCO(ts) + (dt/6)*wCO(ts) + (2/3)*dt*resPenduloCO(2,2) + (1/6)*dt*resPenduloCO(3,2);
        % v{ts + dt}
        wCO(ts + 1) = wCO(ts) + (1/6)*dt*aat(rCO(ts),vCO(ts),titaCO(ts),wCO(ts)) + (2/3)*dt*aat(resResorteCO(2,1), resResorteCO(2,2), resPenduloCO(2,1), resPenduloCO(2,2)) + (1/6)*dt*aat(resResorteCO(3,1), resResorteCO(3,2), resPenduloCO(3,1), resPenduloCO(3,2));

        ts = ts + 1;
    endwhile
    
    
    
    
    
    
    plot(time, titaRK, "r", time, titaTO, "b", time, titaCO, "g") 
    xlabel ("Tiempo t", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    h = legend ("Runge-Kutta", "Tercer orden", "Cuarto orden");
    legend (h, "location", "northeastoutside");
    set (h, "fontsize", 20)
    title ("Péndulo simple", "fontsize", 30);

 
end


