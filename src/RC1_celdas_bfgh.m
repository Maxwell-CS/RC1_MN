%% ============================================================
%  CELDA 2 — b)  ANALISIS GRAFICO   [reemplaza la celda actual]
% ============================================================
% b.1) Hb(Q) y Hs(Q,f(Q)) en [0.001, 0.08] con 500 muestras.
Q_vec  = linspace(0.001, 0.08, 500);
Re_vec = 4*Q_vec/(pi*P.D*P.nu);
f_vec  = nan(size(Q_vec));
conv_v = false(size(Q_vec));
itf_v  = zeros(size(Q_vec));

for i = 1:numel(Q_vec)
    [f_vec(i), itf_v(i), conv_v(i), ~] = friccion_pf(Q_vec(i), P, 1e-10, 100);
end

regimen_ok   = all(Re_vec >= 4000);
convergio_ok = all(conv_v);
fprintf('Regimen turbulento en todo el rango (Re >= 4000): %d  (Re_min = %.0f)\n', ...
        regimen_ok, min(Re_vec));
fprintf('Convergencia del calculo interno de f en las 500 muestras: %d  (iter max = %d)\n', ...
        convergio_ok, max(itf_v));
assert(regimen_ok && convergio_ok, 'Falla el regimen o la convergencia interna.');

Hb_vec = P.H0 + P.a*Q_vec - P.b*Q_vec.^2;
Hs_vec = P.dz + P.C*(f_vec*P.L/P.D + P.KT).*Q_vec.^2;

figure;
plot(Q_vec, Hb_vec, 'b', 'LineWidth', 2); hold on;
plot(Q_vec, Hs_vec, 'r', 'LineWidth', 2);
legend('H_b(Q) bomba','H_s(Q,f(Q)) instalacion','Location','northwest');
title('Altura entregada y altura exigida');
xlabel('Q [m^3/s]'); ylabel('H [m]'); grid on;

% b.2) F(Q) y localizacion de los dos intervalos con cambio de signo.
F_vec = Hb_vec - Hs_vec;

figure;
plot(Q_vec, F_vec, 'g', 'LineWidth', 2); hold on;
yline(0, 'k--');
title('$\mathcal{F}(Q) = H_b(Q) - H_s(Q,f(Q))$', 'Interpreter','latex');
xlabel('Q [m^3/s]'); ylabel('F(Q) [m]'); grid on;

idx = find(F_vec(1:end-1).*F_vec(2:end) < 0);
assert(numel(idx) == 2, 'Se esperaban 2 cambios de signo y se hallaron %d.', numel(idx));
intervalos = [Q_vec(idx)' Q_vec(idx+1)'];        % 2x2 : [l u] por raiz

% Evidencia requerida: valor de F en los extremos de cada intervalo.
evidencia = table( ...
    ["P1";"P1";"P2";"P2"], ...
    [intervalos(1,1); intervalos(1,2); intervalos(2,1); intervalos(2,2)], ...
    [F_vec(idx(1)); F_vec(idx(1)+1); F_vec(idx(2)); F_vec(idx(2)+1)], ...
    'VariableNames', {'Raiz','Q','F_de_Q'})

producto_signos = [F_vec(idx(1))*F_vec(idx(1)+1); F_vec(idx(2))*F_vec(idx(2)+1)];
fprintf('F(l)*F(u) = %+.4e y %+.4e  (ambos < 0: Bolzano justifica biseccion)\n', ...
        producto_signos(1), producto_signos(2));

% Ampliacion cerca de cada cruce.
for r = 1:2
    Qz = linspace(intervalos(r,1), intervalos(r,2), 150);
    Fz = arrayfun(@(q) Fesc_(q,P,1e-10,100), Qz);
    figure;
    plot(Qz, Fz, 'g', 'LineWidth', 2); hold on; yline(0,'k--');
    title(sprintf('Ampliacion cerca de P_%d', r));
    xlabel('Q [m^3/s]'); ylabel('F(Q) [m]'); grid on;
end

% b.3) Denominacion de los equilibrios.
fprintf('\nP1 = equilibrio de MENOR caudal, en [%.6f, %.6f] m^3/s\n', ...
        intervalos(1,1), intervalos(1,2));
fprintf('P2 = equilibrio de MAYOR caudal, en [%.6f, %.6f] m^3/s\n', ...
        intervalos(2,1), intervalos(2,2));
% Ambos satisfacen las ecuaciones, pero eso NO implica que ambos sean
% estables: la estabilidad fisica se decide con el signo de F'(Q_op) en
% el inciso g), no por el hecho de ser raiz.


%% ============================================================
%  CELDA 9 — f)  PROPAGACION DE INCERTIDUMBRES
% ============================================================
% Se resuelve  J(x) dx/dp = -dF/dp  para p en {eps, D}, con las derivadas
% respecto al parametro tomadas a Q y f CONSTANTES.
% D aparece en C (como D^-4), en L/D, en eps/D y en Re: las cuatro
% dependencias estan incluidas en dFdD_.

d_eps = 0.10*P.eps;
d_D   = 0.01*P.D;
raices = sortrows(raices_N);          % [Q f] de P1 y P2, obtenidas en d)

dQde = zeros(2,1);  dQdD = zeros(2,1);
dfde = zeros(2,1);  dfdD = zeros(2,1);
c_eps = zeros(2,1); c_D = zeros(2,1);

for r = 1:2
    Q = raices(r,1);  f = raices(r,2);
    J = Jac_(Q,f,P);
    se = J\(-dFdeps_(Q,f,P));         % dx/d(eps)
    sD = J\(-dFdD_(Q,f,P));           % dx/dD
    dQde(r) = se(1);  dfde(r) = se(2);
    dQdD(r) = sD(1);  dfdD(r) = sD(2);
    c_eps(r) = abs(dQde(r))*d_eps;    % contribucion de eps
    c_D(r)   = abs(dQdD(r))*d_D;      % contribucion de D
end

dQ_lin  = c_eps + c_D;
pct     = 100*dQ_lin./raices(:,1);
dominante = strings(2,1);
for r = 1:2
    if c_D(r) > c_eps(r), dominante(r) = "D"; else, dominante(r) = "eps"; end
end

Tabla_sensibilidad = table(nombres, raices(:,1), dQde, dQdD, c_eps, c_D, dQ_lin, pct, dominante, ...
    'VariableNames', {'Raiz','Q','dQ_deps','dQ_dD','Contrib_eps','Contrib_D', ...
                      'Suma_dQlin','Porc_dQ_Q','Dominante'})

% --- Verificacion por perturbaciones centrales, siguiendo la misma rama ---
pasos = [1e-3 1e-4];
ver = [];
for r = 1:2
    x0 = raices(r,:)';
    for h = pasos
        he = h*P.eps;  hD = h*P.D;

        Pp = P; Pp.eps = P.eps + he;  xp = newton_rc(x0,Pp,1e-12,1e-12,50);
        Pm = P; Pm.eps = P.eps - he;  xm = newton_rc(x0,Pm,1e-12,1e-12,50);
        num_e = (xp(1) - xm(1))/(2*he);

        Pp = P; Pp.D = P.D + hD;      xp = newton_rc(x0,Pp,1e-12,1e-12,50);
        Pm = P; Pm.D = P.D - hD;      xm = newton_rc(x0,Pm,1e-12,1e-12,50);
        num_D = (xp(1) - xm(1))/(2*hD);

        ver = [ver; r h dQde(r) num_e abs(num_e-dQde(r)) ...
                        dQdD(r) num_D abs(num_D-dQdD(r))];   %#ok<AGROW>
    end
end
Tabla_verificacion_sens = array2table(ver, 'VariableNames', ...
    {'Raiz','paso_rel','dQ_deps_analitico','dQ_deps_numerico','dif_eps', ...
     'dQ_dD_analitico','dQ_dD_numerico','dif_D'})

% Interpretacion:
% D domina la sensibilidad en ambas raices, con una contribucion del orden
% de seis veces la de eps. Es consistente con que eps entra solo en A
% (efecto logaritmico y amortiguado), mientras D entra en cuatro lugares y
% el mas fuerte es C ~ D^-4, de modo que un 1% en D mueve las perdidas
% cerca de un 5%. La suma es una cota LINEAL DE PEOR CASO: supone que las
% dos desviaciones ocurren a la vez y en el sentido mas desfavorable. No
% es una desviacion estandar ni un intervalo de confianza, y es distinta
% del error de iteracion, que aqui es del orden de 1e-8 y por lo tanto
% despreciable frente a esta incertidumbre de datos.


%% ============================================================
%  CELDA 10 — g)  INTERPRETACION DE LOS DOS EQUILIBRIOS
% ============================================================
Qg = raices(:,1);  fg = raices(:,2);
Hb_g = zeros(2,1); Hs_g = zeros(2,1); V_g = zeros(2,1); Re_g = zeros(2,1);
cierre = zeros(2,1); Fp_g = zeros(2,1); fp_g = zeros(2,1);
estabilidad = strings(2,1);

for r = 1:2
    Hb_g(r) = Hb_(Qg(r),P);
    Hs_g(r) = Hs_(Qg(r),fg(r),P);
    cierre(r) = abs(Hb_g(r) - Hs_g(r));
    V_g(r)  = 4*Qg(r)/(pi*P.D^2);
    Re_g(r) = Re_(Qg(r),P);
    [Fp_g(r), fp_g(r)] = dF_(Qg(r),P);
    if Fp_g(r) < 0
        estabilidad(r) = "localmente estable";
    elseif Fp_g(r) > 0
        estabilidad(r) = "inestable";
    else
        estabilidad(r) = "la linealizacion no decide";
    end
end

Tabla_equilibrios = table(nombres, Qg, Hb_g, Hs_g, cierre, fg, V_g, Re_g, fp_g, Fp_g, estabilidad, ...
    'VariableNames', {'Punto','Q','Hb','Hs','Cierre_energia','f','V','Re', ...
                      'df_dQ','dF_dQ','Estabilidad'})

fprintf('Cierre de energia: |Hb-Hs| = %.2e y %.2e m  (ambos en el ruido numerico)\n', ...
        cierre(1), cierre(2));
fprintf('Regimen turbulento: Re = %.0f y %.0f, ambos >> 4000\n', Re_g(1), Re_g(2));

% Interpretacion bajo el modelo I dQ/dt = F(Q):
% En P1, F'(Q) > 0. Una perturbacion positiva de caudal hace F > 0, luego
% dQ/dt > 0 y el caudal sigue creciendo, alejandose de P1 hasta alcanzar
% P2. Una perturbacion negativa hace F < 0 y el caudal colapsa hacia el
% reposo. P1 es entonces inestable: es una raiz legitima del sistema, pero
% la instalacion no puede sostenerse en ella.
% En P2, F'(Q) < 0. Una perturbacion positiva produce F < 0, que frena el
% caudal, y una negativa produce F > 0, que lo acelera: en ambos casos el
% sistema regresa a P2. P2 es el punto de operacion localmente estable.
% Notar que la estabilidad de esta iteracion no es lo mismo que la
% estabilidad del metodo de punto fijo del inciso e): alli P2 atraia la
% iteracion por rho < 1, que es una propiedad del algoritmo, no de la
% fisica del sistema.


%% ============================================================
%  CELDA 11 — h)  VARIACION DE LA ALTURA ESTATICA
% ============================================================
% Como F(Q;dz) = Phi(Q) - dz, con Phi(Q) = Hb(Q) - C*B(f(Q))*Q^2
% independiente de dz, todo el barrido consiste en invertir Phi, y
% ademas F'(Q) = Phi'(Q) para cualquier dz.

Qmin = 4000*pi*P.D*P.nu/4;
fprintf('Qmin = %.6e m^3/s  (Re = %.1f)\n', Qmin, Re_(Qmin,P));

dz_vec = 28:0.1:36;
Qg_mesh = linspace(Qmin, 0.08, 4000);
Phi_mesh = arrayfun(@(q) Phi_(q,P), Qg_mesh);

% h.1-h.2) Equilibrios en cada dz y clasificacion de ramas.
res = [];                      % [dz  Q  F'(Q)  estable]
nraices = zeros(size(dz_vec));
for i = 1:numel(dz_vec)
    dz = dz_vec(i);
    G  = Phi_mesh - dz;
    id = find(G(1:end-1).*G(2:end) < 0);
    Qr = [];
    for k = 1:numel(id)
        l = Qg_mesh(id(k));  u = Qg_mesh(id(k)+1);
        for it = 1:200                        % refinamiento por biseccion
            m = (l+u)/2;
            if (Phi_(l,P)-dz)*(Phi_(m,P)-dz) < 0, u = m; else, l = m; end
            if (u-l) <= 1e-14*max(1,abs(m)), break; end
        end
        Qr(end+1) = (l+u)/2;                  %#ok<AGROW>
    end
    if ~isempty(Qr)
        Qr = uniquetol(Qr, 1e-8);             % eliminar duplicados
        for k = 1:numel(Qr)
            Fp = dF_(Qr(k),P);
            res = [res; dz Qr(k) Fp double(Fp<0)];   %#ok<AGROW>
        end
    end
    nraices(i) = numel(Qr);
end

figure; hold on; grid on;
est = res(:,4)==1;
plot(res(est,1),  res(est,2),  '.', 'Color',[0.10 0.45 0.20], 'MarkerSize', 9);
plot(res(~est,1), res(~est,2), '.', 'Color',[0.75 0.15 0.15], 'MarkerSize', 9);
yline(Qmin, 'k--', 'LineWidth', 1.2);
xlabel('\Delta z [m]'); ylabel('Q_{op} [m^3/s]');
title('Ramas de operacion frente a la altura estatica');
legend('rama estable (F'' < 0)','rama inestable (F'' > 0)','Q_{min} (Re = 4000)', ...
       'Location','northeast');

% h.3) Punto de union de las ramas: F'(Qc) = 0 (biseccion sobre F').
lc = Qmin;  uc = 0.08;
for it = 1:200
    mc = (lc+uc)/2;
    if dF_(lc,P)*dF_(mc,P) < 0, uc = mc; else, lc = mc; end
    if (uc-lc) <= 1e-15*max(1,abs(mc)), break; end
end
Qc  = (lc+uc)/2;
dzc = Phi_(Qc,P);
Pc = P; Pc.dz = dzc;
[Fqc, fqc] = Fesc_(Qc, Pc, 1e-12, 200);
fprintf('\nPunto tangente: Qc = %.10f m^3/s , dz_c = %.6f m\n', Qc, dzc);
fprintf('  Comprobacion: F(Qc; dz_c) = %.3e   y   F''(Qc) = %.3e\n', Fqc, dF_(Qc,P));

% h.4) Altura correspondiente a Qmin (salida de una raiz por el dominio).
dz_min = Phi_(Qmin,P);
fprintf('Altura en el limite turbulento: Phi(Qmin) = %.6f m\n', dz_min);

% h.5) Intervalos con 0, 1 o 2 soluciones.
Tabla_conteo = table(dz_vec', nraices', 'VariableNames', {'dz','n_raices'});
for n = 0:2
    v = dz_vec(nraices == n);
    if ~isempty(v)
        fprintf('  %d solucion(es): dz de %.1f a %.1f m (%d valores de la malla)\n', ...
                n, min(v), max(v), numel(v));
    end
end

Tabla_puntos_criticos = table( ...
    ["union tangente de ramas"; "salida por el limite turbulento"], ...
    [Qc; Qmin], [dzc; dz_min], ...
    ["F'(Qc) = 0, las dos raices se funden"; "la raiz de menor caudal cruza Q_min"], ...
    'VariableNames', {'Causa','Q','dz','Descripcion'})

% Las dos causas son distintas: en dz_c las dos raices se funden en un
% cruce TANGENTE, donde F no cambia de signo y por eso una malla de 0.1 m
% no lo localiza (hay que resolver F'(Qc)=0); en dz_min, en cambio, la
% raiz de menor caudal sale del dominio por el borde Re = 4000, sin que
% nada degenere. "Sin raiz turbulenta" no significa "sin flujo posible":
% significa que fuera de ese rango el modelo de Colebrook adoptado ya no
% es aplicable.
