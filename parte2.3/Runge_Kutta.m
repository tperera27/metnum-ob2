function Runge_Kutta (intervals, tmax, dt, r0, v0, tita0, w0)

    %Calculo aceleración
    tita_a = @(tita, v_tita, r, v) ((-2*v*v_tita - 9.81*sin(tita)) / (0.5+r));
    a = @(tita, v_tita, r) ((v_tita**2)*(0.5 + r) + cos(tita)*9.81 - 98.1 * r);
    
    %Phi precalculado
    phi = cell(3,4);
    
    phi{1,1} = 1;
    phi{1,2} = dt / 2;

    phi{2,1} = 1;
    phi{2,2} = -(dt);
    phi{2,3} = 2*dt;

    phi{3,1} = 1;
    phi{3,2} = (dt/6);
    phi{3,3} = (2/3)*dt;
    phi{3,4} = (dt/6);

    %Parámetros
    tao = {1/2, 1, 1};
    alpha = cell(3,3);
    alpha{1,1} = 0;
    alpha{2,1} = 0;
    alpha{2,2} = 0;
    alpha{3,1} = 0;
    alpha{3,2} = 0;
    alpha{3,3} = 0;
    
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
    resResorte = zeros(2,2);
    resPendulo = zeros(2,2);

        
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        % Resorte primer paso
        resResorte(1,1) = phi{1, 1}*r(ts) + phi{1,2}*v(ts);
        resResorte(1,2) = phi{1, 1}*v(ts) + phi{1,2}*a(tita(ts), w(ts), r(ts));

        % Pendulo primer paso
        resPendulo(1,1) = phi{1, 1}*tita(ts) + phi{1,2}*w(ts);
        resPendulo(1,2) = phi{1, 1}*w(ts) + phi{1,2}*tita_a(tita(ts), w(ts), r(ts), v(ts));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        %Resorte segundo paso
        resResorte(2,1) = phi{2,1}*r(ts) + phi{2,2}*v(ts) + phi{2,3}*resResorte(1,2);
        resResorte(2,2) = phi{2,1}*v(ts) + phi{2,2}*a(tita(ts), w(ts), r(ts)) + phi{2,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1));
        
        %Pendulo segundo paso
        resPendulo(2,1) = phi{2,1}*tita(ts) + phi{2,2}*w(ts) + phi{2,3}*resPendulo(1,2);
        resPendulo(2,2) = phi{2,1}*w(ts) + phi{2,2}*tita_a(tita(ts), w(ts), r(ts), v(ts)) + phi{2,3}*tita_a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2));
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Resorte tercer paso
        r(ts + 1) = phi{3,1}*r(ts) + phi{3,2}*v(ts) + phi{3,3}*resResorte(1,2) + phi{3,4}*resResorte(2,2);
        v(ts + 1) = phi{3,1}*v(ts) + phi{3,2}*a(tita(ts), w(ts), r(ts)) + phi{3,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1)) + phi{3,4}*a(resPendulo(2,1), resPendulo(2,2), resResorte(2,1));

        % Pendulo tercer paso
        tita(ts + 1) = phi{3,1}*tita(ts) + phi{3,2}*w(ts) + phi{3,3}*resPendulo(1,2) + phi{3,4}*resPendulo(2,2);
        w(ts + 1) = phi{3,1}*w(ts) + phi{3,2}*tita_a(tita(ts), w(ts), r(ts), v(ts)) + phi{3,3}*tita_a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2)) + phi{3,4}*tita_a(resPendulo(2,1), resPendulo(2,2), resResorte(2,1), resResorte(2,2));

        ts = ts + 1;
    endwhile
    
    plot(time, tita, "r");
    xlabel ("Tiempo", "fontsize", 20);
    ylabel ("Posición angular", "fontsize", 20);
    title ("Péndulo resorte - Método de Runge-Kutta de tercer orden", "fontsize", 30);
end


