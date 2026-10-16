%% TP3 - Identification des parametres geometriques d'un robot 3R plan
clear; clc; close all; rng(1);

p_true = [0.2 0.3 0.2 0.1 -0.1 0.08 0.1 0.03 0.025 0.6]';
names  = {'l0','l1','l2','l3','o1','o2','o3','xc','yc','ac'};
idx    = [2 3 4 6 7 8 9 10];      % parametres identifiables (l0 et o1 fixes)

%% 2. Verification du jacobien (differences finies)
th  = 2*pi*rand(3,1) - pi;
dp  = 1e-6*randn(10,1);
dX  = modele_geom(th, p_true+dp) - modele_geom(th, p_true);
err = norm(dX - jacobien(th,p_true)*dp);
fprintf('Verif jacobien : ||dX - J.dp|| = %.2e (||dX|| = %.2e)\n\n', err, norm(dX));

%% 3. Jacobien 30x10, rang
N  = 10;
Q  = 2*pi*rand(N,3) - pi;
[~, J] = empile(Q, p_true);
fprintf('Jacobien %dx%d : rang = %d\n', size(J,1), size(J,2), rank(J));
fprintf('Valeurs singulieres :\n'); disp(svd(J)');
Jr = J(:, idx);
fprintf('Avec l0 et o1 imposes : rang = %d sur %d colonnes, cond = %.2e\n\n', ...
        rank(Jr), numel(idx), cond(Jr));

%% 4. Identification sans erreur de mesure
Xmes = empile(Q, p_true);
p0   = p_true;  p0(idx) = p_true(idx) + 0.1*randn(numel(idx),1);   % estimation initiale
[p_id, cost] = identifier(Q, Xmes, p0, idx, 100);
fprintf('--- Identification sans bruit ---\n');
disp(table(names', p_true, p0, p_id, 'VariableNames', {'param','vrai','initial','identifie'}));
fprintf('Erreur max : %.2e  (cout final %.2e)\n\n', max(abs(p_id-p_true)), cost(end));

%% 5. Robustesse
sig_q = 1e-3;   % bruit articulaire (rad)
sig_x = 1e-3;   % bruit cartesien (m)  (et rad pour alpha)
nMC   = 100;

% --- 5a. Influence du nombre de mesures
Nlist = [4 5 6 8 10 15 20 30 50 100];
err_a = zeros(size(Nlist));
for k = 1:numel(Nlist)
    e = zeros(nMC,1);
    for m = 1:nMC
        Qk = 2*pi*rand(Nlist(k),3) - pi;
        [Xn, Qn] = mesures_bruitees(Qk, p_true, sig_q, sig_x);
        p0 = p_true; p0(idx) = p_true(idx) + 0.05*randn(numel(idx),1);
        pid = identifier(Qn, Xn, p0, idx, 50);
        e(m) = norm(pid(idx) - p_true(idx));
    end
    err_a(k) = mean(e);
end
figure; loglog(Nlist, err_a, 'o-','LineWidth',1.5); grid on;
xlabel('Nombre de mesures'); ylabel('Erreur moyenne ||p_{id}-p||');
title('5a : influence du nombre de mesures');

% --- 5b. Influence du choix des mesures
cas = {'theta1 seul','theta1 et theta2','theta2 et theta3','theta1,2,3'};
msk = {[1 0 0],[1 1 0],[0 1 1],[1 1 1]};
Nb  = 30;
for c = 1:numel(cas)
    e = zeros(nMC,1);
    Qk0 = (2*pi*rand(Nb,3) - pi) .* repmat(msk{c}, Nb, 1);
    Jc  = empile(Qk0, p_true);  [~,Jc] = empile(Qk0, p_true);
    for m = 1:nMC
        Qk = (2*pi*rand(Nb,3) - pi) .* repmat(msk{c}, Nb, 1);
        [Xn, Qn] = mesures_bruitees(Qk, p_true, sig_q, sig_x);
        p0 = p_true; p0(idx) = p_true(idx) + 0.05*randn(numel(idx),1);
        pid = identifier(Qn, Xn, p0, idx, 50);
        e(m) = norm(pid(idx) - p_true(idx));
    end
    fprintf('5b : %-18s rang(J)=%d  cond=%.1e  erreur moyenne=%.2e\n', ...
        cas{c}, rank(Jc(:,idx)), cond(Jc(:,idx)), mean(e));
end

% --- 5c. Influence de l'estimation initiale (minima locaux)
nT = 200;  Nc = 30;
Qk = 2*pi*rand(Nc,3) - pi;
[Xn, Qn] = mesures_bruitees(Qk, p_true, sig_q, sig_x);
ok = 0; fin = zeros(nT,1); errs = zeros(nT,1);
for t = 1:nT
    p0 = p_true;
    p0(idx) = [0.5*rand(3,1); 2*pi*rand(2,1)-pi; 0.2*randn(2,1); 2*pi*rand-pi];
    [pid, cost] = identifier(Qn, Xn, p0, idx, 200);
    fin(t) = cost(end);  errs(t) = norm(pid(idx)-p_true(idx));
    ok = ok + (errs(t) < 0.05);
end
fprintf('5c : %d/%d tirages convergent vers le minimum global\n', ok, nT);
figure; semilogy(sort(fin),'.'); grid on;
xlabel('tirages tries'); ylabel('cout final'); title('5c : minima locaux');

%% ================= Fonctions locales =================
function [Y, J] = empile(Q, p)
% Empile les modeles (3N x 1) et jacobiens (3N x 10) pour N mesures
N = size(Q,1);
Y = zeros(3*N,1);  J = zeros(3*N,10);
for i = 1:N
    Y(3*i-2:3*i)    = modele_geom(Q(i,:), p);
    J(3*i-2:3*i,:)  = jacobien(Q(i,:), p);
end
end

function [Xn, Qn] = mesures_bruitees(Q, p, sq, sx)
% Mesures simulees avec bruit gaussien (articulaire + cartesien)
N  = size(Q,1);
Qn = Q + sq*randn(N,3);            % angles mesures bruites
X  = empile(Q, p);                 % mesures cartesiennes exactes (angles vrais)
Xn = X + sx*randn(3*N,1);
end

function [p, cost] = identifier(Q, Xmes, p0, idx, nIter)
% Moindres carres non lineaires (Gauss-Newton + pseudo-inverse)
p = p0(:);  cost = zeros(nIter,1);
for k = 1:nIter
    [Y, J] = empile(Q, p);
    e = Xmes - Y;
    e(3:3:end) = atan2(sin(e(3:3:end)), cos(e(3:3:end)));   % residu angulaire modulo 2pi
    cost(k) = 0.5*(e'*e);
    dp = pinv(J(:,idx)) * e;
    p(idx) = p(idx) + dp;
    if norm(dp) < 1e-12, cost = cost(1:k); break; end
end
end
