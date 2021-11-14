

function Tercer_Orden(intervals, tmax, dt, r0, v0, tita0, w0)
    
    phi = cell(3,4);
    
    phi{1,1} = 1;
    phi{1,2} = dt / 3;

    phi{2,1} = 1;
    phi{2,2} = 0;
    phi{2,3} = (2*dt) / 3;

    phi{3,1} = 1;
    phi{3,2} = (dt/4);
    phi{3,3} = 0;
    phi{3,4} = ((3*dt)/4);

    tao = {1/3, 2/3, 1};
    alpha = zeros(3,3);
    alpha(1,1) = 1;
    alpha(2,1) = -2/3;
    alpha(2,2) = 2/3;
    alpha(3,1) = 1/3;
    alpha(3,2) = -2/3;
    alpha(3,3) = 1/3;
    
    
    time = 0:dt:tmax;
    
    
    r = zeros(1, intervals);
    r(1) = r0;
    
    
    v = zeros(1, intervals);
    v(1) = v0;
    
    tita = zeros(1, intervals);
    tita(1) = tita0;
    
    w = zeros(1, intervals);
    w(1) = w0;
    
    resResorte = zeros(2,2);
    resPendulo = zeros(2,2);
    
    % Aceleración del resorte (radial)
    ar = @(t) (((w(t))^2)*(0.5 + r(t)) + (cos(tita(t))*9.81) - (98.1 * r(t)));
    art = @(rt,vt,titat,wt) (((wt^2)*(0.5 + rt) + (cos(titat))*9.81) - (98.1 * rt));
    
    % Aceleración del pendulo (angular)
    aa = @(t) (-( (2*r(t)*tita(t)) + (9.81*sin(tita(t))) ) / (0.5 + r(t)) );
    aat = @(rt,vt,titat,wt) (-( (2*rt*titat) + (9.81*sin(titat)) ) / (0.5 + rt) );
   
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 4 y 5
        
        % Resorte primer paso
        % u{ts + dt/3}
        resResorte(1,1) = phi{1,1}*r(ts) + phi{1,2}*v(ts) + (1/2)*((dt/3)^2)*(alpha(1,1))*ar(ts);
        % v{ts + dt/3}
        resResorte(1,2) = phi{1,1}*v(ts) + phi{1,2}*ar(ts);
        
        % Pendulo primer paso
        % u{ts + dt/3}
        resPendulo(1,1) = phi{1,1}*tita(ts) + phi{1,2}*w(ts) + (1/2)*((dt/3)^2)*(alpha(1,1))*aa(ts);
        % v{ts + dt/3}
        resPendulo(1,2) = phi{1,1}*w(ts) + phi{1,2}*aa(ts);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 8 y 9
        
        % Resorte segundo paso
        % u{ts + 2dt/3}
        resResorte(2,1) = phi{2,1}*r(ts) + phi{2,2}*v(ts) + phi{2,3}*resResorte(1,2) + (0.5)*((2*dt/3)^2)*( (alpha(2,1)*ar(ts)) + (alpha(2,2)*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) );
        % v{ts + 2dt/3}
        resResorte(2,2) = phi{2,1}*v(ts) + phi{2,2}*ar(ts) + phi{2,3}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2));

        % Pendulo segundo paso
        % u{ts + 2dt/3}
        resPendulo(2,1) = phi{2,1}*tita(ts) + phi{2,2}*w(ts) + phi{2,3}*resPendulo(1,2) + (0.5)*((2*dt/3)^2)*( (alpha(2,1)*aa(ts)) + (alpha(2,2)*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) );
        % v{ts + 2dt/3}
        resPendulo(2,2) = phi{2,1}*w(ts) + phi{2,2}*aa(ts) + phi{2,3}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 13 y 14
        
        % Resorte segundo paso
        % u{ts + dt}
        r(ts + 1) = phi{3,1}*r(ts) + phi{3,2}*v(ts) + phi{3,3}*resResorte(1,2) + phi{3,4}*resResorte(2,2) + (0.5)*(dt^2)*( (alpha(3,1)*ar(ts)) + (alpha(3,2)*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) + (alpha(3,3)*art(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2))) );
        % v{ts + dt}
        v(ts + 1) = phi{3,1}*v(ts) + phi{3,2}*ar(ts) + phi{3,3}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2)) + phi{3,4}*art(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2));
        
        % Pendulo segundo paso
        % u{ts + dt}
        tita(ts + 1) = phi{3,1}*tita(ts) + phi{3,2}*w(ts) + phi{3,3}*resPendulo(1,2) + phi{3,4}*resPendulo(2,2) + (0.5)*(dt^2)*( (alpha(3,1)*aa(ts)) + (alpha(3,2)*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) + (alpha(3,3)*aat(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2))) );
        % v{ts + dt}
        w(ts + 1) = phi{3,1}*w(ts) + phi{3,2}*aa(ts) + phi{3,3}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2)) + phi{3,4}*aat(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2));
        
        ts = ts + 1;
    endwhile
        
    plot(time, tita, "b");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo resorte - Método explicito de tercer orden", "fontsize", 30);
end