%% ============================================================
%  CELDA FINAL — FUNCIONES
%  Va como ULTIMA celda de codigo del MLX (lo exige la guia base:
%  "adjuntelas integramente al final del MLX").
%  Contiene:
%    (A) Funciones del laboratorio de clase, transcritas SIN cambios.
%    (B) Rutinas del caso, construidas sobre los algoritmos de clase
%        y adaptadas a los criterios (7)-(8) del enunciado.
% ============================================================


%% ---------- (A) FUNCIONES DE CLASE (Laboratorio S4) ----------
% Transcritas literalmente del material del curso. Se adjuntan por
% exigencia de la guia base. Observacion importante: ambas rutinas de
% Newton forman inv(J) de manera explicita y usan una norma infinito
% mezclada como criterio de parada. El enunciado de este RC pide lo
% contrario en los dos puntos (resolver J*s=-F sin invertir, y error
% relativo POR COMPONENTE con escalas Qs y fs), por lo que en el
% inciso d) se emplea la version adaptada newton_rc de la seccion (B).

function z = newton2v(F,x0,Tol)
    syms x y
    dF = jacobian(F(x,y),[x,y]);
    dF_fun = matlabFunction(dF);
    dr = 1;
    i = 0;
    z = [i x0' dr];
    while dr > Tol
        x1 = x0 - inv(dF_fun(x0(1),x0(2)))*F(x0(1),x0(2));
        dr = norm(x1 - x0,'inf')/norm(x1,'inf');
        x0 = x1;
        i = i + 1;
        z = [z; i x1' dr];
    end
end

function z = newtonVec(F,x0,Maxiter)
    n = length(x0);
    syms x [n 1]
    J_sym = jacobian(F(x),x);
    J = matlabFunction(J_sym, 'Vars',{x});
    error = NaN;
    z = [0 x0' error];
    for i = 1:Maxiter
        x1 = x0 - inv(J(x0))*F(x0);
        error = norm(x1 - x0,'inf')/norm(x1,'inf');
        x0 = x1;
        z = [z; i x1' error];
    end
end

function z = pfijos(G,x0,Tol)
    dr = 1;
    i = 0;
    z = [i x0' dr];
    while dr > Tol
        x1 = G(x0);
        dr = norm(x1-x0,'inf')/norm(x1,'inf');
        x0 = x1;
        i = i + 1;
        z = [z; i x1' dr];
    end
end


%% ---------- (B) RUTINAS DEL CASO ----------

% ----- Bloque basico del modelo -----

function r = Re_(Q,P)
    r = 4*Q./(pi*P.D*P.nu);
end

function A = A_(Q,f,P)
    A = P.eps/(3.7*P.D) + 2.51./(Re_(Q,P).*sqrt(f));
end

function B = B_(f,P)
    B = f*P.L/P.D + P.KT;
end

function v = N_(Q,P)
    v = P.H0 + P.a*Q - P.b*Q.^2 - P.dz;
end

function v = Hb_(Q,P)
    v = P.H0 + P.a*Q - P.b*Q.^2;
end

function v = Hs_(Q,f,P)
    v = P.dz + P.C*B_(f,P).*Q.^2;
end

function v = F1_(Q,f,P)
    v = N_(Q,P) - P.C*B_(f,P).*Q.^2;
end

function v = F2_(Q,f,P)
    v = 1./sqrt(f) + 2*log10(A_(Q,f,P));
end

function v = R_(Q,f,P)
    % Residuo adimensional, Ec. (8)
    v = max(abs(F1_(Q,f,P))/P.Href, abs(F2_(Q,f,P)));
end

function Jm = Jac_(Q,f,P)
    % Jacobiano analitico J = d(F1,F2)/d(Q,f), con la otra variable fija.
    t  = 2.51/(Re_(Q,P)*sqrt(f));
    A  = A_(Q,f,P);
    Jm = [ P.a - 2*P.b*Q - 2*P.C*B_(f,P)*Q ,  -P.C*(P.L/P.D)*Q^2 ; ...
           -2*t/(log(10)*Q*A)              ,  -0.5*f^(-3/2) - t/(log(10)*f*A) ];
end

function [ok,causa] = dominio_(Q,f,P)
    % Dominio real declarado: Q>0, f>0, Re>=4000, 0<A<1.
    ok = false;  causa = "-";
    if ~isfinite(Q) || ~isfinite(f), causa = "valor no finito";     return; end
    if ~(Q > 0),                     causa = "Q <= 0";              return; end
    if ~(f > 0),                     causa = "f <= 0";              return; end
    if Re_(Q,P) < 4000,              causa = "Re < 4000";           return; end
    A = A_(Q,f,P);
    if ~(A > 0 && A < 1),            causa = "A fuera de (0,1)";    return; end
    ok = true;
end

% ----- Iteracion interna: factor de friccion por punto fijo -----

function [f,it,ok,causa] = friccion_pf(Q,P,tol,maxit)
    % Resuelve Colebrook (Ec. 3) para Q fijo, por punto fijo
    % f_{k+1} = [-2 log10 A(Q,f_k)]^(-2), semilla 0.02.
    % Declara convergencia solo si el error relativo de f Y |F2| <= tol.
    f = 0.02;  it = 0;  ok = false;  causa = "-";
    if ~(Q > 0),          causa = "Q <= 0";     f = NaN;  return; end
    if Re_(Q,P) < 4000,   causa = "Re < 4000";  f = NaN;  return; end
    for k = 1:maxit
        A = A_(Q,f,P);
        if ~(A > 0 && A < 1)
            causa = "A fuera de (0,1)";  f = NaN;  it = k;  return;
        end
        fn = 1/(-2*log10(A))^2;
        er = abs(fn - f)/max(abs(fn), 1e-12*P.fs);
        f  = fn;  it = k;
        if er <= tol && abs(F2_(Q,f,P)) <= tol
            ok = true;  return;
        end
    end
    causa = "no converge en f";
end

function [Fq,f,itf,ok,causa] = Fesc_(Q,P,tol,maxit)
    % Funcion escalar F(Q) = F1(Q, f(Q)) de la Ec. (6).
    [f,itf,ok,causa] = friccion_pf(Q,P,tol,maxit);
    if ~ok,  Fq = NaN;  return;  end
    Fq = F1_(Q,f,P);
end

function v = Phi_(Q,P)
    % Altura estatica que hace F(Q)=0:  Phi(Q) = Hb(Q) - C*B(f(Q))*Q^2.
    % Como F(Q;dz) = Phi(Q) - dz, el barrido del inciso h) es invertir Phi.
    [f,~,ok,~] = friccion_pf(Q,P,1e-12,200);
    if ~ok,  v = NaN;  return;  end
    v = Hb_(Q,P) - P.C*B_(f,P)*Q^2;
end

function [Fp,fp] = dF_(Q,P)
    % Ec. (13):  f'(Q) = -F2_Q/F2_f ,  F'(Q) = F1_Q + F1_f*f'(Q).
    [f,~,ok,~] = friccion_pf(Q,P,1e-12,200);
    if ~ok,  Fp = NaN;  fp = NaN;  return;  end
    J  = Jac_(Q,f,P);
    fp = -J(2,1)/J(2,2);
    Fp = J(1,1) + J(1,2)*fp;
end

% ----- Derivadas respecto a los parametros (inciso f) -----

function d = dFdeps_(Q,f,P)
    % dF/d(eps) a Q,f constantes.  eps solo entra en A.
    A = A_(Q,f,P);
    d = [ 0 ; 2/(log(10)*A) * (1/(3.7*P.D)) ];
end

function d = dFdD_(Q,f,P)
    % dF/dD a Q,f constantes.  D entra en C (D^-4), en L/D, en eps/D y en Re.
    %   C*B = 8*f*L/(g*pi^2*D^5) + 8*KT/(g*pi^2*D^4)
    %   d(C*B)/dD = -40*f*L/(g*pi^2*D^6) - 32*KT/(g*pi^2*D^5)
    % y en F2, el termino 2.51/(Re*sqrt(f)) = 2.51*pi*D*nu/(4*Q*sqrt(f)) es
    % proporcional a D, de modo que d(t)/dD = t/D.
    A  = A_(Q,f,P);
    t  = 2.51/(Re_(Q,P)*sqrt(f));
    dF1 = Q^2*( 40*f*P.L/(P.g*pi^2*P.D^6) + 32*P.KT/(P.g*pi^2*P.D^5) );
    dF2 = 2/(log(10)*A) * ( -P.eps/(3.7*P.D^2) + t/P.D );
    d = [ dF1 ; dF2 ];
end

% ----- Biseccion (algoritmo de clase S3, criterios del enunciado) -----

function [raiz,T,salida,nEval] = biseccion_rc(l,u,P,tolE,tolR,N,verb)
    % Algoritmo de biseccion de S3 (Bolzano + division a la mitad),
    % con los criterios del enunciado:  E_k=(u-l)/(2|m|)<=tolE  y
    % |F(m)|/Href<=tolR.  La iteracion interna de f usa tol 1e-10.
    % Si una evaluacion de F falla, NO se usa para actualizar el intervalo.
    if nargin < 7, verb = true; end
    tolInt = 1e-10;  maxInt = 100;
    nEval  = 0;  raiz = NaN;  salida = "limite de iteraciones";
    datos  = zeros(0,9);

    [Fl,~,itl,okl] = Fesc_(l,P,tolInt,maxInt);  nEval = nEval + 1;
    [Fu,~,itu,oku] = Fesc_(u,P,tolInt,maxInt);  nEval = nEval + 1;
    if ~(okl && oku)
        salida = "dominio invalido en un extremo";
        T = array2table(datos,'VariableNames', ...
            {'k','l','u','m','f_m','F_m','E_k','R_m','IterInternas'});
        return
    end
    if Fl*Fu >= 0
        salida = "sin cambio de signo (Bolzano no se cumple)";
        T = array2table(datos,'VariableNames', ...
            {'k','l','u','m','f_m','F_m','E_k','R_m','IterInternas'});
        return
    end
    if verb
        fprintf('  Bolzano: F(%.6f) = %+.6e , F(%.6f) = %+.6e  ->  producto %+.3e < 0\n', ...
                l, Fl, u, Fu, Fl*Fu);
    end

    for k = 1:N
        m = (l+u)/2;
        [Fm,fm,itf,okm] = Fesc_(m,P,tolInt,maxInt);  nEval = nEval + 1;
        if ~okm
            salida = "dominio invalido en m";  break
        end
        Ek = (u-l)/(2*abs(m));
        Rm = abs(Fm)/P.Href;
        datos(end+1,:) = [k l u m fm Fm Ek Rm itf];   %#ok<AGROW>
        raiz = m;
        if Fm == 0
            salida = "raiz exacta";  break
        end
        if Ek <= tolE && Rm <= tolR
            salida = "convergio";  break
        end
        if Fl*Fm < 0
            u = m;
        else
            l = m;  Fl = Fm;
        end
    end
    T = array2table(datos,'VariableNames', ...
        {'k','l','u','m','f_m','F_m','E_k','R_m','IterInternas'});
end

% ----- Newton para el sistema acoplado (adaptado del laboratorio S4) -----

function [x,T,salida,nF,nJ] = newton_rc(x0,P,tole,tolR,Maxiter)
    % Adaptacion de newtonVec/newton2v del laboratorio S4 a lo que pide
    % este RC.  Tres cambios deliberados respecto de la version de clase:
    %   1) el paso se obtiene resolviendo  J(x_k) s_k = -F(x_k)  con el
    %      operador \ , sin formar inv(J);
    %   2) el Jacobiano es el analitico Jac_, no el simbolico;
    %   3) el criterio de parada es el (7)-(8) por componente, con las
    %      escalas Qs y fs, en lugar de la norma infinito mezclada.
    Q = x0(1);  f = x0(2);
    nF = 0;  nJ = 0;
    datos = zeros(0,5);
    salida = "limite de iteraciones";

    [ok,causa] = dominio_(Q,f,P);
    if ~ok
        salida = "dominio invalido en la semilla: " + causa;
        x = [Q;f];
        T = array2table(datos,'VariableNames',{'k','Q_k','f_k','e_k','R_k'});
        return
    end
    datos(end+1,:) = [0 Q f NaN R_(Q,f,P)];  nF = nF + 1;

    for k = 1:Maxiter
        Fv = [F1_(Q,f,P); F2_(Q,f,P)];  nF = nF + 1;
        J  = Jac_(Q,f,P);               nJ = nJ + 1;
        if rcond(J) < eps
            salida = "Jacobiano mal condicionado";  break
        end
        s  = J\(-Fv);                   % sin formar inv(J)
        Qn = Q + s(1);
        fn = f + s(2);

        [ok,causa] = dominio_(Qn,fn,P);
        if ~ok
            salida = "dominio invalido: " + causa;
            datos(end+1,:) = [k Qn fn NaN NaN];   %#ok<AGROW>
            break
        end

        ek = max( abs(Qn-Q)/max(abs(Qn),1e-12*P.Qs), ...
                  abs(fn-f)/max(abs(fn),1e-12*P.fs) );
        Q = Qn;  f = fn;
        Rk = R_(Q,f,P);  nF = nF + 1;
        datos(end+1,:) = [k Q f ek Rk];   %#ok<AGROW>

        if ek <= tole && Rk <= tolR
            salida = "convergio";  break
        end
    end
    x = [Q;f];
    T = array2table(datos,'VariableNames',{'k','Q_k','f_k','e_k','R_k'});
end
