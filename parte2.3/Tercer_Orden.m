function Tercer_Orden(intervals, tmax, dt, r0, v0, tita0, w0)
    
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
    alphaTO = zeros(3,3);
    alphaTO(1,1) = 1;
    alphaTO(2,1) = -2/3;
    alphaTO(2,2) = 2/3;
    alphaTO(3,1) = 1/3;
    alphaTO(3,2) = -2/3;
    alphaTO(3,3) = 1/3;
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
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
        
    plot(time, titaTO, "b");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo resorte - Método explicito de tercer orden", "fontsize", 30);
end