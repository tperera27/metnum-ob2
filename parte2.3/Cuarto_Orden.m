function Cuarto_Orden(intervals, tmax, dt, r0, v0, tita0, w0)
    
    %Phi precalculado
    phi = cell(4,5);
    
    phi{1,1} = 1;
    phi{1,2} = dt / 3;

    phi{2,1} = 1;
    phi{2,2} = dt/8;
    phi{2,3} = (3/8)*dt;

    phi{3,1} = 1;
    phi{3,2} = (dt/2);
    phi{3,3} = -(3/2)*dt;
    phi{3,4} = 2*dt;
    
    phi{4,1} = 1;
    phi{4,2} = (dt/6);
    phi{4,3} = 0;
    phi{4,4} = (2/3)*dt;
    phi{4,5} = (dt/6);
    
    %Parámetros
    tao = {1/3, 1/2, 1};
    alpha = cell(4,4);
    alpha{1,1} = 1;
    alpha{2,1} = -3/5;
    alpha{2,2} = 3/5;
    alpha{3,1} = 3/5;
    alpha{3,2} = 3/5;
    alpha{3,3} = 0;
    alpha{4,1} = 0;
    alpha{4,2} = 0;
    alpha{4,3} = 0;
    alpha{4,4} = 0;
    
    %Vector con los tiempos a graficar
    time = 0:dt:tmax;
    
    %Posición radial
    r = zeros(1, intervals);
    r(1) = r0;
    
    %Velocidad radial
    v = zeros(1, intervals);
    v(1) = v0;
    
    %Posición angular
    tita = zeros(1, intervals);
    tita(1) = tita0;
    
    %Velocidad angular
    w = zeros(1, intervals);
    w(1) = w0;
    
    %Resultados intermedios
    resResorte = zeros(3,2);
    resPendulo = zeros(3,2);
    
    %Calculo aceleración del resorte (radial)
    ar = @(t) (((w(t))^2)*(0.5 + r(t)) + (cos(tita(t))*9.81) - (98.1 * r(t)));
    art = @(rt,vt,titat,wt) (((wt^2)*(0.5 + rt) + (cos(titat))*9.81) - (98.1 * rt));
    
    %Calculo aceleración del pendulo (angular)
    aa = @(t) (-( (2*r(t)*tita(t)) + (9.81*sin(tita(t))) ) / (0.5 + r(t)) );
    aat = @(rt,vt,titat,wt) (-( (2*rt*titat) + (9.81*sin(titat)) ) / (0.5 + rt) );
   
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 4 y 5
        
        % Resorte primer paso
        % u{ts + dt/3}
        resResorte(1,1) = phi{1,1}*r(ts) + phi{1,2}*v(ts) + (1/2)*((dt/3)^2)*alpha{1,1}*ar(ts);
        % v{ts + dt/3}
        resResorte(1,2) = phi{1,1}*v(ts) + phi{1,2}*ar(ts);
        
        % Pendulo primer paso
        % u{ts + dt/3}
        resPendulo(1,1) = phi{1,1}*tita(ts) + phi{1,2}*w(ts) + (1/2)*((dt/3)^2)*(alpha{1,1})*aa(ts);
        % v{ts + dt/3}
        resPendulo(1,2) = phi{1,1}*w(ts) + phi{1,2}*aa(ts);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuaciones 8 y 9
        
        % Resorte segundo paso
        % u{ts + dt/2}
        resResorte(2,1) = phi{2,1}*r(ts) + phi{2,2}*v(ts) + phi{2,3}*resResorte(1,2) + (0.5)*((dt/2)^2)*( (alpha{2,1}*ar(ts)) + (alpha{2,2}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) );
        % v{ts + dt/2}
        resResorte(2,2) = phi{2,1}*v(ts) + phi{2,2}*ar(ts) + phi{2,3}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2));

        % Pendulo segundo paso
        % u{ts + dt/2}
        resPendulo(2,1) = phi{2,1}*tita(ts) + phi{2,2}*w(ts) + phi{2,3}*resPendulo(1,2) + (0.5)*((dt/2)^2)*( (alpha{2,1}*aa(ts)) + (alpha{2,2}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) );
        % v{ts + dt/2}
        resPendulo(2,2) = phi{2,1}*w(ts) + phi{2,2}*aa(ts) + phi{2,3}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuacion 28 y 29
        
        % Resorte tercer paso
        % u{ts + dt}
        resResorte(3,1) = phi{3,1}*r(ts) + phi{3,2}*v(ts) + phi{3,3}*resResorte(1,2) + phi{3,4}*resResorte(2,2) + (0.5)*(dt^2)*( (alpha{3,1}*ar(ts)) + (alpha{3,2}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) + (alpha{3,3}*art(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2))) );
        % v{ts + dt}
        resResorte(3,2) = phi{3,1}*v(ts) + phi{3,2}*ar(ts) + phi{3,3}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2)) + phi{3,4}*art(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2));

        % Pendulo tercer paso
        % u{ts + dt}
        resPendulo(3,1) = phi{3,1}*tita(ts) + phi{3,2}*w(ts) + phi{3,3}*resPendulo(1,2) + phi{3,4}*resPendulo(2,2) + (0.5)*(dt^2)*( (alpha{3,1}*aa(ts)) + (alpha{3,2}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2))) + (alpha{3,3}*aat(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2))) );
        % v{ts + dt}
        resPendulo(3,2) = phi{3,1}*w(ts) + phi{3,2}*aa(ts) + phi{3,3}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2)) + phi{3,4}*aat(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Cuarto paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ecuacion 30 y 31
        
        % Resorte cuarto paso
        % u{ts + dt}
        r(ts + 1) = phi{4,1}*r(ts) + phi{4,2}*v(ts) + phi{4,3}*resResorte(1,2) + phi{4,4}*resResorte(2,2) + phi{4,5}*resResorte(3,2);
        % v{ts + dt}
        v(ts + 1) = phi{4,1}*v(ts) + phi{4,2}*ar(ts) + phi{4,3}*art(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2)) + phi{4,4}*art(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2)) + phi{4,5}*art(resResorte(3,1), resResorte(3,2), resPendulo(3,1), resPendulo(3,2));

        % Pendulo cuarto paso
        % u{ts + dt}
        tita(ts + 1) = phi{4,1}*tita(ts) + phi{4,2}*w(ts) + phi{4,3}*resPendulo(1,2) + phi{4,4}*resPendulo(2,2) + phi{4,5}*resPendulo(3,2);
        % v{ts + dt}
        w(ts + 1) = phi{4,1}*w(ts) + phi{4,2}*aa(ts) + phi{4,3}*aat(resResorte(1,1), resResorte(1,2), resPendulo(1,1), resPendulo(1,2)) + phi{4,4}*aat(resResorte(2,1), resResorte(2,2), resPendulo(2,1), resPendulo(2,2)) + phi{4,5}*aat(resResorte(3,1), resResorte(3,2), resPendulo(3,1), resPendulo(3,2));

        ts = ts + 1;
    endwhile
    
    disp(tita)
    plot(time, tita, "g");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo resorte - Método explicito de cuarto orden", "fontsize", 30);
end