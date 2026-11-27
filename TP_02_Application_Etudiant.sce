// Executer depuis le dossier contenant les deux fichiers Scilab.
clear; clc;
// Le script retrouve ses fichiers meme si le dossier courant est different.
chdir(get_absolute_file_path("TP_02_Application_Etudiant.sce"));
exec("TP_02_Outils_EDO.sci",-1);

function z=f_test(t,y)
    z=1;
endfunction
[t,y,c]=integrer_edo(f_test,0,1,2,10,"euler");
if max(abs(y-(2+t)))>1d-12 then error("Test constant echoue"); end
mprintf("Test constant Euler : OK, cout=%d\n",c);
// ===== Controles Q1 / Q2 : Euler, Heun, RK4 =====
methodes_test=["euler" "heun" "rk4"];
couts_test=[1 2 4]; // evaluations de f par pas

// --- Q1 : un pas pour y'=-2y, y0=1, h=0,1 ---
function z=f_lin(t,y)
z=-2*y;
endfunction
attendu_q1=[0.8 0.82 0.8187333333333333];

for i=1:3
methode=methodes_test(i);
[t,y,c]=integrer_edo(f_lin,0,0.1,1,1,methode);

if abs(y(2)-attendu_q1(i))>1d-12 then
error("Controle Q1 echoue : "+methode);
end
if c<>couts_test(i) then
error("Cout Q1 incorrect : "+methode);
end

mprintf("Q1 %s : y1=%.12f, cout=%d : OK\n",methode,y(2),c);
end

// --- Q2 : y'=1, y(0)=2 sur [0,1], N=10 (f_test est defini plus haut) ---
for i=1:3
methode=methodes_test(i);
[t,y,c]=integrer_edo(f_test,0,1,2,10,methode);

if max(abs(y-(2+t)))>1d-12 then
error("Test constant echoue : "+methode);
end
if size(t,"*")<>11 | size(y,"*")<>11 then
error("Nombre de valeurs incorrect : "+methode);
end
if abs(t(1))>1d-12 | abs(t($)-1)>1d-12 then
error("Maillage incorrect : "+methode);
end
if c<>10*couts_test(i) then
error("Cout incorrect : "+methode);
end

mprintf("Q2 constant %s : OK, cout=%d\n",methode,c);
end

// --- Q2 : controle du temps intermediaire, y'=t, y(0)=0, un pas h=0,2 ---
function z=f_temps(t,y)
z=t;
endfunction
attendu_temps=[0 0.02 0.02];

for i=1:3
methode=methodes_test(i);
[t,y,c]=integrer_edo(f_temps,0,0.2,0,1,methode);

if abs(y($)-attendu_temps(i))>1d-12 then
error("Test temps echoue : "+methode);
end

mprintf("Q2 temps %s : y1=%.4f : OK\n",methode,y($));
end
// ===== Fin des controles Q1 / Q2 =====
function z=f_decay(t,y)
    z=-15*y;
endfunction
Ns=[5 10 20 40 80];
methodes=["euler"]; // TODO : ajouter "heun" et "rk4" une fois completes.
resultats=[]; // methode_id : 1=euler, 2=heun, 3=rk4 (ordre fixe).
// Conserver cet ordre en ajoutant les methodes.
csvWrite([1;2;3],"TP2_methodes.csv",";",".","%d", ..
    "methode_id (1=euler, 2=heun, 3=rk4)");
for m=1:size(methodes,"*")
    E=zeros(1,size(Ns,"*")); cout=E;
    for j=1:size(Ns,"*")
        [t,y,c]=integrer_edo(f_decay,0,1,1,Ns(j),methodes(m));
        // TODO A : remplacer %nan par l'erreur maximale aux noeuds.
        E(j)=%nan;
        cout(j)=c;
    end
    for j=1:size(Ns,"*")
        p=%nan;
        if j<size(Ns,"*") then p=ordre_edo(E(j),E(j+1)); end
        resultats=[resultats; m Ns(j) 1/Ns(j) E(j) cout(j) p];
    end
end
csvWrite(resultats,"TP2_convergence.csv",";",".","%.12e", ..
    "methode_id;N;h;Emax;nfev;p_vers_ligne_suivante");
disp(resultats);
// Des NaN dans Emax signalent TODO A non termine.
// TODO B : tableau a budget egal (Euler 400, Heun 200, RK4 100 pas).
// TODO C : trajectoires et erreurs ; exclure NaN/zeros du trace log.
// CSV final : ajouter diagnostic stable/instable et ordre interpretable.
// Derniere ligne de chaque methode : p=NaN est normal (pas de suivant).

// Application thermique : un modele pedagogique, pas des mesures.
Cth=1000; H=10; P=100; Ta=20; theta0=20;
function z=f_therm(t,theta)
    z=(P-H*(theta-Ta))/Cth;
endfunction
function theta=theta_ref(t,Hval)
    eq=Ta+P/Hval;
    theta=eq+(theta0-eq)*exp(-Hval*t/Cth);
endfunction
// TODO D : horizon 1000 s, N=100,10,4 ; comparer les trois methodes.
// TODO E : reference H=8,10,12 W/K a temps fixes ; interpreter les ecarts.
// Conseils graphiques dans l'interface Scilab :
// plot2d(h,E,logflag="ll") pour un tableau d'erreurs positives.
// Ajouter titres, legendes, unites et identifier les calculs instables.
