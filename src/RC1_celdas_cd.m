%% ============================================================
%  CELDA 1 — a)  SISTEMA NO LINEAL   [reemplaza la celda actual]
% ============================================================
clearvars; clc; close all;

% --- Datos del caso (Tabla del enunciado, unidades SI) ---
P.H0  = 30;        % m        Altura a caudal nulo
P.a   = 300;       % s/m^2    Coeficiente lineal de la bomba
P.b   = 3750;      % s^2/m^5  Coeficiente cuadratico de la bomba
P.dz  = 32;        % m        Desnivel entre superficies libres
P.D   = 0.20;      % m        Diametro interior
P.L   = 200;       % m        Longitud de la tuberia
P.eps = 4.5e-5;    % m        Rugosidad absoluta
P.KT  = 5;         % -        Suma de coeficientes de perdidas locales
P.nu  = 1.0e-6;    % m^2/s    Viscosidad cinematica
P.g   = 9.81;      % m/s^2    Aceleracion de la gravedad

% --- Agrupacion del enunciado ---
P.C = 8/(P.g*pi^2*P.D^4);     % s^2/m^5

% --- Escalas y referencias de los criterios (7)-(8) ---
P.Qs   = 0.01;     % m^3/s
P.fs   = 0.02;     % -
P.Href = 30;       % m

% --- Tolerancias y limites declarados (reproducibilidad) ---
tol_e   = 1e-6;    % error relativo por componente, Ec. (7)
tol_R   = 1e-8;    % residuo adimensional, Ec. (8)
tol_int = 1e-10;   % iteracion interna de f (Colebrook)
maxit_int = 100;
maxit_bis = 100;
maxit_new = 50;

% --- Verificacion de unidades ---
% F1 esta en metros:   [C*B(f)*Q^2] = (s^2/m^5)(-)(m^3/s)^2 = m
% F2 es adimensional:  [f^-1/2] = (-)  y  [A] = (-) porque
%     [eps/(3.7 D)] = m/m = -   y   [2.51/(Re sqrt(f))] = -
fprintf('Verificacion de unidades\n');
fprintf('  [C]      = s^2/m^5   -> C = %.4f\n', P.C);
fprintf('  [C*B*Q^2]= m         -> con Q=0.05, f=0.02: %.4f m\n', ...
        P.C*(0.02*P.L/P.D + P.KT)*0.05^2);
fprintf('  [A]      = adimensional (rugosidad relativa + termino de Re)\n\n');

% Nota: F = (F1,F2)^T depende de DOS variables (Q,f); la funcion escalar
% F(Q) = F1(Q,f(Q)) depende de una sola, porque f se obtiene resolviendo
% F2(Q,f)=0 internamente para cada Q.


%% ============================================================
%  CELDA 3 — c)  BISECCION CON ITERACION INTERNA   [reemplaza la actual]
% ============================================================
% Se aplica biseccion en LOS DOS intervalos hallados en b).
% La iteracion interna de f usa punto fijo con semilla 0.02 y tol 1e-10;
% si falla, esa evaluacion de F no se usa para actualizar el intervalo.

assert(exist('intervalos','var')==1, ...
    'Ejecute antes la celda b): debe existir la matriz "intervalos" (2x2).');

nombres  = ["P1"; "P2"];
raiz_bis = zeros(2,1);
fbis     = zeros(2,1);
iter_bis = zeros(2,1);
evals_bis= zeros(2,1);
estado_bis = strings(2,1);
Tbis     = cell(2,1);

for r = 1:2
    l = intervalos(r,1);  u = intervalos(r,2);
    fprintf('\n--- Biseccion en %s, intervalo [%.6f, %.6f] ---\n', nombres(r), l, u);
    [m,T,salida,nEv] = biseccion_rc(l,u,P,1e-6,1e-8,maxit_bis);
    Tbis{r}       = T;
    raiz_bis(r)   = m;
    iter_bis(r)   = height(T);
    evals_bis(r)  = nEv;
    estado_bis(r) = salida;
    if ~isempty(T), fbis(r) = T.f_m(end); else, fbis(r) = NaN; end
    fprintf('  Estado: %s | iteraciones: %d | evaluaciones de F: %d\n', ...
            salida, iter_bis(r), nEv);
    fprintf('  Raiz: Q = %.10f m^3/s , f = %.10f\n', raiz_bis(r), fbis(r));
    disp(T)
end

Tabla_biseccion = table(nombres, raiz_bis, fbis, iter_bis, evals_bis, estado_bis, ...
    'VariableNames', {'Raiz','Q','f','Iteraciones','EvalF','Estado'})

% Por que la tolerancia interna es mas estricta:
% el valor f(Q) del punto fijo entra en F(Q), que es la funcion sobre la
% que se aplica biseccion. Un error en f se propaga a F y puede alterar su
% SIGNO cerca de la raiz, provocando una actualizacion incorrecta del
% intervalo. Con tol interna 1e-10 frente a 1e-8 del residuo externo, el
% error heredado de f queda dos ordenes por debajo del criterio de parada,
% de modo que el signo de F(m) es confiable en todas las iteraciones.


%% ============================================================
%  CELDA 4 — d)  NEWTON PARA EL SISTEMA ACOPLADO   [reemplaza la actual]
% ============================================================
% Jacobiano analitico (derivadas parciales, con la otra variable FIJA):
%
%   dF1/dQ = a - 2 b Q - 2 C B(f) Q
%   dF1/df = -C (L/D) Q^2
%   dF2/dQ = -2 t / (ln10 * Q * A)                 con t = 2.51/(Re sqrt(f))
%   dF2/df = -0.5 f^(-3/2) - t / (ln10 * f * A)
%
% Acoplamiento: F1 necesita f porque las perdidas dependen de la friccion,
% y F2 necesita Q porque el numero de Reynolds depende del caudal. Ninguna
% de las dos se puede resolver aislada.

% --- Verificacion del Jacobiano con derivacion simbolica ---
syms Qv fv real
Cv  = 8/(P.g*pi^2*P.D^4);
Bv  = fv*P.L/P.D + P.KT;
Rev = 4*Qv/(pi*P.D*P.nu);
Av  = P.eps/(3.7*P.D) + 2.51/(Rev*sqrt(fv));
F1v = P.H0 + P.a*Qv - P.b*Qv^2 - P.dz - Cv*Bv*Qv^2;
F2v = 1/sqrt(fv) + 2*log10(Av);
J_sym = jacobian([F1v; F2v], [Qv, fv]);
J_num = matlabFunction(J_sym, 'Vars', {Qv, fv});
dif_J = norm(J_num(0.02,0.025) - Jac_(0.02,0.025,P), inf);
fprintf('Diferencia entre el Jacobiano analitico y el simbolico: %.3e\n', dif_J);
assert(dif_J < 1e-9, 'El Jacobiano analitico no coincide con el simbolico.');

% --- Newton desde las dos semillas del enunciado ---
semillas_N = [0.005 0.030; 0.060 0.020];
raices_N   = zeros(2,2);
iter_N     = zeros(2,1);
nF_N       = zeros(2,1);
nJ_N       = zeros(2,1);
estado_N   = strings(2,1);
TN         = cell(2,1);

for s = 1:2
    [x,T,salida,nF,nJ] = newton_rc(semillas_N(s,:)', P, tol_e, tol_R, maxit_new);
    TN{s} = T;  raices_N(s,:) = x';  estado_N(s) = salida;
    iter_N(s) = height(T)-1;  nF_N(s) = nF;  nJ_N(s) = nJ;
    fprintf('\nSemilla (Q0,f0) = (%.3f, %.3f)  ->  %s\n', ...
            semillas_N(s,1), semillas_N(s,2), salida);
    fprintf('  Q = %.10f m^3/s , f = %.10f , iteraciones externas = %d\n', ...
            x(1), x(2), iter_N(s));
    disp(T)
end

Tabla_newton = table(semillas_N(:,1), semillas_N(:,2), raices_N(:,1), raices_N(:,2), ...
    iter_N, nF_N, nJ_N, estado_N, ...
    'VariableNames', {'Q0','f0','Q','f','IterExternas','EvalF','EvalJ','Estado'})

% --- Comparacion con biseccion ---
raices_N_ord = sortrows(raices_N);
dif_metodos  = abs(raices_N_ord(:,1) - sort(raiz_bis));
Tabla_comparacion_cd = table(nombres, sort(raiz_bis), raices_N_ord(:,1), dif_metodos, ...
    'VariableNames', {'Raiz','Q_biseccion','Q_newton','Diferencia'})

% Por que las iteraciones externas no miden el costo total:
% cada iteracion de Newton evalua F, evalua J y resuelve un sistema 2x2,
% pero trabaja con (Q,f) simultaneos. Cada evaluacion de F en biseccion
% arrastra una iteracion interna completa de punto fijo sobre Colebrook
% (del orden de 7-10 pasos). Por eso hay que comparar EVALUACIONES DE
% FUNCION, no iteraciones externas: ver la tabla de sintesis final.
