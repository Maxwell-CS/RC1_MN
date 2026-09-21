%% ============================================================
%  CELDA 12 — SINTESIS Y VERIFICACION
%  (va despues de h) y antes de la celda de funciones)
% ============================================================
% Reune los tres metodos. Una salida por dominio invalido NO se reporta
% como convergencia: aparece con su propio estado.

metodo = strings(0,1); semilla = strings(0,1); raizAlc = strings(0,1);
Qfinal = []; iters = []; evalsF = []; Rfinal = []; salida = strings(0,1);

% --- Biseccion (inciso c) ---
for r = 1:2
    metodo(end+1,1)  = "Biseccion";
    semilla(end+1,1) = sprintf("[%.5f, %.5f]", intervalos(r,1), intervalos(r,2));
    raizAlc(end+1,1) = nombres(r);
    Qfinal(end+1,1)  = raiz_bis(r);
    iters(end+1,1)   = iter_bis(r);
    evalsF(end+1,1)  = evals_bis(r);
    Rfinal(end+1,1)  = Tbis{r}.R_m(end);
    salida(end+1,1)  = estado_bis(r);
end

% --- Newton (inciso d) ---
for s = 1:2
    metodo(end+1,1)  = "Newton";
    semilla(end+1,1) = sprintf("(%.3f, %.3f)", semillas_N(s,1), semillas_N(s,2));
    [~, ir] = min(abs(raices_N(s,1) - raices(:,1)));
    if estado_N(s) == "convergio", raizAlc(end+1,1) = nombres(ir);
    else,                          raizAlc(end+1,1) = "ninguna"; end
    Qfinal(end+1,1) = raices_N(s,1);
    iters(end+1,1)  = iter_N(s);
    evalsF(end+1,1) = nF_N(s);
    Rfinal(end+1,1) = R_(raices_N(s,1), raices_N(s,2), P);
    salida(end+1,1) = estado_N(s);
end

% --- Punto fijo (inciso e; usa las variables ya calculadas alli) ---
for s = 1:numel(estado_pf)
    metodo(end+1,1)  = "Punto fijo";
    semilla(end+1,1) = sprintf("(%.3f, %.3f)", semillas_pf(s,1), semillas_pf(s,2));
    raizAlc(end+1,1) = destino_pf(s);
    Qfinal(end+1,1)  = Qfin_pf(s);
    iters(end+1,1)   = iter_pf(s);
    evalsF(end+1,1)  = 2*iter_pf(s);      % g1 y g2 por iteracion
    Rfinal(end+1,1)  = Rfin_pf(s);
    salida(end+1,1)  = estado_pf(s);
end

Tabla_comparativa_metodos = table(metodo, semilla, raizAlc, Qfinal, iters, evalsF, Rfinal, salida, ...
    'VariableNames', {'Metodo','Semilla_o_Intervalo','Raiz_alcanzada','Q_final', ...
                      'Iteraciones','Eval_funcion','Residuo_final','Estado_salida'})

fprintf('\n--- Conclusiones ---\n');
fprintf('Equilibrio localmente estable: %s, con Q = %.6f m^3/s (F'' = %+.2f < 0).\n', ...
        nombres(Fp_g<0), Qg(Fp_g<0), Fp_g(Fp_g<0));
fprintf('Equilibrio inestable: %s, con Q = %.6f m^3/s (F'' = %+.2f > 0).\n', ...
        nombres(Fp_g>0), Qg(Fp_g>0), Fp_g(Fp_g>0));
fprintf('Parametro dominante en la sensibilidad: D, con %.1f%% y %.1f%% de la\n', ...
        100*c_D(1)/dQ_lin(1), 100*c_D(2)/dQ_lin(2));
fprintf('  desviacion lineal total en P1 y P2 respectivamente.\n');

% Conclusion:
% Bajo el modelo dinamico adoptado, el unico punto de operacion sostenible
% es P2 (mayor caudal), donde F'(Q) < 0. P1 satisface las ecuaciones pero
% es inestable, y ningun metodo iterativo lo "descarta" por eso: Newton lo
% alcanza sin dificultad desde una semilla cercana, mientras que el punto
% fijo no llega a el porque rho(JG) > 1. Son dos nociones de estabilidad
% distintas y no deben confundirse.
% El parametro dominante en la sensibilidad es el diametro D, que entra en
% cuatro dependencias del modelo; una tolerancia de fabricacion del 1% en D
% produce mas incertidumbre en el caudal que una variacion del 10% en la
% rugosidad.
% Para seleccionar una bomba real faltaria informacion que este modelo no
% contiene: margen de NPSH disponible frente al requerido (cavitacion),
% curva de eficiencia y potencia al eje en el punto de operacion, limites
% de caudal minimo y maximo continuo del fabricante, envejecimiento de la
% tuberia (eps crece con los anos y desplaza el punto de operacion), y el
% comportamiento en transitorios, que el modelo cuasiestacionario adoptado
% no describe.
