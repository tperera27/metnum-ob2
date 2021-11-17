function Runge_Kutta (intervals, tmax, dt, r0, v0, tita0, w0)

    %Calculo aceleración
    tita_a = @(tita, w, r, v) ((-2*v*w - 9.81*sin(tita)) / (0.5+r));
    a = @(tita, w, r) ((w^2)*(0.5 + r) + cos(tita)*9.81 - 98.1 * r);
    
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
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
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
    
    plot(time, titaRK, "r");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo resorte - Método de Runge-Kutta de tercer orden", "fontsize", 30);
end


