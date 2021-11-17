function Cuarto_Orden(intervals, tmax, dt, r0, v0, tita0, w0)
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
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
    
    plot(time, titaCO, "g");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo resorte - Método explicito de cuarto orden", "fontsize", 30);
end