function Runge_Kutta (intervals, tmax, dt, x0, v0, tita0, w0)
    phi = cell(3,4);
    tao1 = 1/2;
    tao2 = 1;


    tita_a = @(tita, v_tita, r, v) ((-2*v*v_tita - 9.81*sin(tita)) / (0.5+r))
    a = @(tita, v_tita, r) ((v_tita**2)*(0.5 + r) + cos(tita)*9.81 - 98.1 * r)

    phi{1,1} = 1;
    phi{1,2} = dt / 2;

    phi{2,1} = 1;
    phi{2,2} = -(dt);
    phi{2,3} = 2*dt;

    phi{3,1} = 1;
    phi{3,2} = (dt/6);
    phi{3,3} = (2/3)*dt;
    phi{3,4} = (dt/6);

    tao = {1/2, 1, 1};
    
    time = 0:dt:tmax;
    x = zeros(1, intervals);
    x(1) = x0;
    v = zeros(1, intervals);
    v(1) = v0;
    tita = zeros(1, intervals);
    tita(1) = tita0;
    w = zeros(1, intervals);
    w(1) = w0;
    resResorte = zeros(2,2);
    resPendulo = zeros(2,2);


    alpha = cell(3,3);
    alpha{1,1} = 0;
    alpha{2,1} = 0;
    alpha{2,2} = 0;
    alpha{3,1} = 0;
    alpha{3,2} = 0;
    alpha{3,3} = 0;
        
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            % Resorte primer paso
            resResorte(1,1) = phi{1, 1}*x(ts) + phi{1,2}*v(ts);
            resResorte(1,2) = phi{1, 1}*v(ts) + phi{1,2}*a(tita(ts), w(ts), x(ts));

            % Pendulo primer paso
            resPendulo(1,1) = phi{1, 1}*tita(ts) + phi{1,2}*w(ts);
            resPendulo(1,2) = phi{1, 1}*w(ts) + phi{1,2}*tita_a(tita(ts), w(ts), x(ts), v(ts));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Resorte segundo paso
            resResorte(2,1) = phi{2,1}*x(ts) + phi{2,2}*v(ts) + phi{2,3}*resResorte(1,2);
            resResorte(2,2) = phi{2,1}*v(ts) + phi{2,2}*a(tita(ts), w(ts), x(ts)) + phi{2,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1));
            
            %Pendulo segundo paso
            resPendulo(2,1) = phi{2,1}*tita(ts) + phi{2,2}*w(ts) + phi{2,3}*resPendulo(1,2);
            resPendulo(2,2) = phi{2,1}*w(ts) + phi{2,2}*tita_a(tita(ts), w(ts), x(ts), v(ts)) + phi{2,3}*tita_a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2));
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Resorte tercer paso
            x(ts + 1) = phi{3,1}*x(ts) + phi{3,2}*v(ts) + phi{3,3}*resResorte(1,2) + phi{3,4}*resResorte(2,2);
            v(ts + 1) = phi{3,1}*v(ts) + phi{3,2}*a(tita(ts), w(ts), x(ts)) + phi{3,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1)) + phi{3,4}*a(resPendulo(2,1), resPendulo(2,2), resResorte(2,1));

            % Pendulo tercer paso
            tita(ts + 1) = phi{3,1}*tita(ts) + phi{3,2}*w(ts) + phi{3,3}*resPendulo(1,2) + phi{3,4}*resPendulo(2,2);
            w(ts + 1) = phi{3,1}*w(ts) + phi{3,2}*tita_a(tita(ts), w(ts), x(ts), v(ts)) + phi{3,3}*tita_a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2)) + phi{3,4}*tita_a(resPendulo(2,1), resPendulo(2,2), resResorte(2,1), resResorte(2,2));

        ts = ts + 1;
    endwhile
    
    plot(time, tita,"r");
    xlabel ("t");
    ylabel ("u(t)");
    title ("");
end


