

function TercerOrden(intervals, tmax, dt, x0, v0, tita0, w0)
    phi = cell(3,4);
    tao1 = 1/3;
    tao2 = 2/3;


    a = @(tita, v_tita, r) ((v_tita**2)*(0.5 + r) + cos(tita)*9.81 - 98.1 * r)
    

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

    alpha = zeros(3,3);
    alpha(1,1) = 1;
    alpha(2,1) = -2/3;
    alpha(2,2) = 2/3;
    alpha(3,1) = 1/3;
    alpha(3,2) = -2/3;
    alpha(3,3) = 1/3;

    
    ts = 1;
    while (ts < intervals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Primer Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Resorte primer paso
        resResorte(1,1) = phi{1, 1}*x(ts) + phi{1,2}*v(ts) + calculadorAlphasResorte(x, tita, w, v,ts, dt, tao{1},alpha(1,1),alpha(1,2),alpha(1,3), resPendulo, resResorte);
        resResorte(1,2) = phi{1, 1}*v(ts) + phi{1,2}*a(tita(ts), w(ts), x(ts));
        
        % Pendulo primer paso
        resPendulo(1,1) = phi{1, 1}*tita(ts) + phi{1,2}*w(ts) + calculadorAlphasPendulo(x, tita, w, v,ts, dt, tao{1},0,0,0,resPendulo ,resResorte);
        resPendulo(1,2) = phi{1, 1}*w(ts) + phi{1,2}*calculadorAceleracionPendulo(tita(ts), w(ts), x(ts), v(ts));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Segundo Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Resorte segundo paso
        resResorte(2,1) = phi{2,1}*x(ts) + phi{2,2}*v(ts) + phi{2,3}*resResorte(1,2) + calculadorAlphasResorte(x, tita, w, v,ts, dt, tao{2},alpha(2,1),alpha(2,2),alpha(2,3), resPendulo, resResorte);
        resResorte(2,2) = phi{2,1}*v(ts) + phi{2,2}*a(tita(ts), w(ts), x(ts)) + phi{2,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2));

        % Pendulo segundo paso
        resPendulo(2,1) = phi{2,1}*tita(ts) + phi{2,2}*w(ts) + phi{2,3}*resResorte(1,2) + calculadorAlphasPendulo(x, tita, w, v,ts, dt, tao{2}, alpha(2,1), alpha(2,2),alpha(2,3), resPendulo, resResorte);
        resPendulo(2,2) = phi{2,1}*w(ts) + phi{2,2}*calculadorAceleracionPendulo(tita(ts), w(ts), x(ts), v(ts)) + phi{2,3}*calculadorAceleracionPendulo(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Tercer Paso        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Resorte tercer paso
        x(ts + 1) = phi{3,1}*x(ts) + phi{3,2}*v(ts) + phi{3,3}*resResorte(1,2) + phi{3,4}*resResorte(2,2) + calculadorAlphasResorte(x, tita, w, v,ts, dt, tao{3},alpha(3,1),alpha(3,2),alpha(3,3), resPendulo, resResorte);
        v(ts + 1) = phi{3,1}*v(ts) + phi{3,2}*a(tita(ts), w(ts), x(ts)) + phi{3,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1)) + phi{3,4}*a(resPendulo(2,1), resPendulo(2,2), resResorte(2,1));
        
        tita(ts + 1) = phi{3,1}*tita(ts) + phi{3,2}*w(ts) + phi{3,3}*resPendulo(1,2) + phi{3,4}*resPendulo(2,2) + calculadorAlphasPendulo(x, tita, w, v,ts, dt, tao{3}, alpha(3,1), alpha(3,2), alpha(3,3), resPendulo, resResorte);
        w(ts + 1) = phi{3,1}*w(ts) + phi{3,2}*a(tita(ts), w(ts), x(ts), v(ts)) + phi{3,3}*a(resPendulo(1,1), resPendulo(1,2), resResorte(1,1), resResorte(1,2)) + phi{3,4}*a(resPendulo(2,1), resPendulo(2,2), resResorte(2,1), resResorte(2,2));
        
        ts = ts + 1;
    endwhile
        
    plot(time, tita, "b");
    xlabel ("t");
    ylabel ("u(t)");
    title ("");
end

function res = calculadorAceleracionResorte(tita, v_tita, r) 
    res = ((v_tita**2)*(0.5 + r) + cos(tita)*9.81 - 98.1 * r)
    return;
end

function res = calculadorAceleracionPendulo(tita, v_tita, r, v) 
    res = ((-2*v*v_tita - 9.81*sin(tita)) / (0.5+r))
    return;
end

function res1 = calculadorAlphasResorte(x, tita, w, v, ts, dt, tao, alpha1=1, alpha2=0, alpha3=0, Pendulo, Resorte) 
    res1 = ( 0.5*((tao*dt)**2)*( alpha1*calculadorAceleracionResorte(tita(ts), w(ts), x(ts)) + alpha2*calculadorAceleracionResorte(Pendulo(1,1), Pendulo(1,2), Resorte(1,1)) + alpha3*calculadorAceleracionResorte(Pendulo(2,1), Pendulo(2,2), Resorte(2,1)) ));
    return;
end
function res2 = calculadorAlphasPendulo(x, tita, w, v, ts, dt, tao, alpha1=1, alpha2=0, alpha3=0, Pendulo, Resorte) 
    res2 = ( 0.5*((tao*dt)**2)*( (alpha1*calculadorAceleracionPendulo(tita(ts), w(ts), x(ts), v(ts))) + (alpha2*calculadorAceleracionPendulo(Pendulo(1,1), Pendulo(1,2), Resorte(1,1), Resorte(1,2))) +(alpha3*calculadorAceleracionPendulo(Pendulo(2,1), Pendulo(2,2), Resorte(2,1), Resorte(2,2)))))
    return;
end