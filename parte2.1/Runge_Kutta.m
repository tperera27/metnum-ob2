function Runge_Kutta (dt, x0, v0)
    phi = cell(3,4);
    tao1 = 1/2;
    tao2 = 1;


    a = @(u) (-(100)*(1+10*u*u)*u);
    

    phi{1,1} = @(t) (1);
    phi{1,2} = @(t) (t);

    phi{2,1} = @(t) (1);
    phi{2,2} = @(t) (-(t*(t + (tao1 * dt))) / (3 * tao1 * dt));
    phi{2,3} = @(t) ((t*t) / (tao1*dt));

    phi{3,1} = @(t) (1);
    phi{3,2} = @(t) ((t*(2*t*t - 3*tao1*dt*t -3*tao2*dt*t + 6*tao1*tao2*dt*dt)) / (6*tao1*tao2*dt*dt));
    phi{3,3} = @(t) (-(t*t*(-2*t + 3*tao2*dt) / (6*tao1*(-tao2+tao1)*dt*dt)));
    phi{3,4} = @(t) ((t*t*(-2*t + 3*tao1*dt) / (6*tao2*(-tao2+tao1)*dt*dt)));

    disp(phi{1,1}(0))

    tao = {1/2, 1, 1};
    resultados = zeros(3,2);

    alpha = cell(3,3);
    alpha{1,1} = 0;
    alpha{2,1} = 0;
    alpha{2,2} = 0;
    alpha{3,1} = 0;
    alpha{3,2} = 0;
    alpha{3,3} = 0;

    disp(phi{1, 1}(tao{1}*dt));

    for i=1:3
        resultados(i,1) = phi{i, 1}(tao{i}*dt)*x0 + phi{i,2}(tao{i}*dt)*v0;
        resultados(i,2) = phi{i, 1}(tao{i}*dt)*v0 + phi{i,2}(tao{i}*dt)*a(x0);

        if (i == 2)
            resultados(i,1) = resultados(i, 1) + phi{i, 3}(dt*tao{i}) * resultados(1,2) + ((tao{i}*dt * tao{i}*dt) / 2) * (alpha{i,1}*a(x0) + alpha{i,i}*a(resultados(i-1,1)));

            resultados(i,2) = resultados(i,2) + phi{i,3}(tao{i}*dt) * a(resultados(i-1, 1));
        endif
        if (i == 3)
            resultados(i,1) = resultados(i, 1) + phi{i, 3}(dt*tao{i}) * resultados(1,2)+ phi{i, 4}(dt*tao{i}) * resultados(2,2) + ((tao{i}*dt * tao{i}*dt) / 2) * (alpha{i,i-2}*a(x0) + alpha{i,i-1}*a(resultados(i-2,1)) + alpha{i,i}*a(resultados(i-1,1)));

            resultados(i,2) = resultados(i,2) + phi{i,3}(tao{i}*dt) * a(resultados(i-2, 1)) + phi{i,4}(tao{i}*dt) * a(resultados(i-1, 1));
        endif
    endfor
    disp(resultados);

    ploting = cell(2,30);
    for i=1:30
        ploting{1,i} = i;
        ploting{2,i} = i*5;
    endfor

    plot(ploting{1, : }, ploting{2, :});
    
end


