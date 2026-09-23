%% APP2
% BOIF1302
% DESM1210

clc
clear
close all

%% Paramètres
Isource      = "lenna.bmp";   % "cman.tif", "irm.tif", "mandrill.tif", "crest.bmp"
taille_cible = 256;

nbits_QS   = 5;  
nbits_DPCM = 5;  
taille_bloc = 2; 
nbits_VQ   = 5;

Tableau_Delta = [];
Tableau_PSNR = [];

%% Chargement
[Isource, IsourceMod] = charger_image(Isource);

%% Redimensionnement
Iredim_PPV = redim_PPV(IsourceMod, taille_cible);
Iredim_BiL = redim_bilineaire(IsourceMod, taille_cible);

%% Compression
% QS
for Delta_QS = 0:0.01:1
    [IDecoder_QS,   PSNR_QS]   = compression_QS(Iredim_PPV, nbits_QS, Delta_QS);
    Tableau_Delta = [Tableau_Delta, Delta_QS];
    Tableau_PSNR = [Tableau_PSNR, PSNR_QS];
end
[Value, idx_QS] = max(Tableau_PSNR(:));
[IDecoder_QS,   PSNR_QS]   = compression_QS(Iredim_PPV, nbits_QS, Tableau_Delta(idx_QS));

% DPCM
Tableau_Delta = [];
Tableau_PSNR = [];
for Delta_DPCM = 0:0.01:1
    [IDecoder_DPCM, PSNR_DPCM] = compression_DPCM(Iredim_BiL, nbits_DPCM, Delta_DPCM);
    Tableau_Delta = [Tableau_Delta, Delta_DPCM];
    Tableau_PSNR = [Tableau_PSNR, PSNR_DPCM];
end
[Value, idx_DPCM] = max(Tableau_PSNR(:));
[IDecoder_DPCM, PSNR_DPCM] = compression_DPCM(Iredim_BiL, nbits_DPCM, Tableau_Delta(idx_DPCM));

% Vectoriel
[IDecoder_VQ,   PSNR_VQ]   = compression_VQ(Iredim_PPV, taille_bloc, nbits_VQ);

%% Affichage
afficher_image(1, Isource,       "ISource "    + taille_str(Isource))
afficher_image(2, IsourceMod,    "ISourceMod " + taille_str(IsourceMod))
afficher_image(3, Iredim_PPV,    "IRedim PPV " + taille_str(Iredim_PPV))
afficher_image(4, Iredim_BiL,    "IRedim BiL " + taille_str(Iredim_BiL))
afficher_image(5, IDecoder_QS,   "IDecoder QS "         + 2^nbits_QS   + " niveaux, PSNR = " + PSNR_QS)
afficher_image(6, IDecoder_DPCM, "IDecoder DPCM "       + 2^nbits_DPCM + " niveaux, PSNR = " + PSNR_DPCM)
afficher_image(7, IDecoder_VQ,   "IDecoder Vectoriel "  + 2^nbits_VQ   + " codewords, PSNR = " + PSNR_VQ)


%% ======================================================================
%%  FONCTIONS
%% ======================================================================

%% --- Chargement --------------------------------------------------------
function [Isource, Igris] = charger_image(Isource)
    Isource = imread(Isource);
    if size(Isource, 3) == 3
        Igris = mean(double(Isource), 3);   % (R+G+B)/3
    else
        Igris = double(Isource);
    end
end

%% --- Redimensionnement plus proche voisin ------------------------------
function Iredim = redim_PPV(I, taille_cible)
    [L1, C1] = size(I);
    if ~(L1 > taille_cible && C1 > taille_cible)
        Iredim = I;
        return
    end
    L2 = taille_cible;  C2 = taille_cible;
    Iredim = zeros(L2, C2);
    for l = 1 : L2
        for c = 1 : C2
            ll = max(floor(l*L1/L2), 1);
            cc = max(floor(c*C1/C2), 1);
            Iredim(l, c) = I(ll, cc);
        end
    end
end

%% --- Redimensionnement bilinéaire --------------------------------------
function Iredim = redim_bilineaire(I, taille_cible)
    [L1, C1] = size(I);
    if ~(L1 > taille_cible && C1 > taille_cible)
        Iredim = I;
        return
    end
    L2 = taille_cible;  C2 = taille_cible;
    Iredim = zeros(L2, C2);
    for l = 1 : L2
        for c = 1 : C2
            ll = min(max(floor(l*L1/L2), 1), L1-1);
            cc = min(max(floor(c*C1/C2), 1), C1-1);
            dl = l*L1/L2 - ll;
            dc = c*C1/C2 - cc;
            X = (1-dc)*I(ll,   cc) + dc*I(ll,   cc+1);
            Y = (1-dc)*I(ll+1, cc) + dc*I(ll+1, cc+1);
            Iredim(l, c) = floor((1-dl)*X + dl*Y);
        end
    end
end

%% --- Quantificateur scalaire -------------------------------------------
function [Idec, PSNR] = compression_QS(I, nbits, Delta)
    Q = genere_quantif(nbits, Delta);
    m     = mean(I(:));
    sigma = std(I(:));

    Inorm  = (I - m) / sigma;
    Icode  = Q93(Inorm, Q);
    InormQ = Q93_1(Icode, Q);
    Idec   = InormQ*sigma + m;

    PSNR = calcul_PSNR(I, Idec);
end

%% --- DPCM --------------------------------------------------------------
function [Idec, PSNR] = compression_DPCM(I, nbits, Delta)
    Q = genere_quantif(nbits, Delta);
    [L, C] = size(I);

    % Prédiction et erreur de prédiction
    xpred = zeros(L, C);
    for l = 1 : L
        for c = 1 : C
            xpred(l, c) = predire_pixel(I, l, c);
        end
    end
    erreur = I - xpred;

    % Normalisation, quantification, reconstruction
    m     = mean(erreur(:));
    sigma = std(erreur(:));
    eq    = Q93_1(Q93((erreur - m)/sigma, Q), Q);
    Idec  = (eq*sigma + m) + xpred;

    PSNR = calcul_PSNR(I, Idec);
end

function xpred = predire_pixel(I, l, c)
    if l == 1 && c == 1
        xpred = 128;
    elseif l == 1
        xpred = I(l, c-1);
    elseif c == 1
        xpred = I(l-1, c);
    else
        p1 = 0.5*I(l-1, c)   + 0.5*I(l, c-1);
        p2 = 0.5*I(l-1, c-1) + 0.5*I(l, c-1);
        p3 = 0.5*I(l-1, c-1) + 0.5*I(l-1, c);
        xpred = median([p1, p2, p3]);
    end
end

%% --- Quantification vectorielle ----------------------------------------
function [Idec, PSNR] = compression_VQ(I, taille_bloc, nbits)
    [L, C] = size(I);
    M = 2^nbits;

    blocs        = decoupe_blocs(I, taille_bloc);        % 1. découpage
    dictionnaire = LBG_splitting(blocs, M);              % 2. dictionnaire
    index        = encode_VQ(blocs, dictionnaire);       % 3. encodage
    blocs_dec    = dictionnaire(:, index);               % 4. décodage
    Idec = recombine_blocs(blocs_dec, taille_bloc, L, C);% 5. reconstruction

    PSNR = calcul_PSNR(I, Idec);
end

function blocs = decoupe_blocs(I, taille_bloc)
    [L, C] = size(I);
    nb_blocs = (L/taille_bloc) * (C/taille_bloc);
    blocs = zeros(taille_bloc^2, nb_blocs);
    idx = 1;
    for l = 1 : taille_bloc : L
        for c = 1 : taille_bloc : C
            b = I(l:l+taille_bloc-1, c:c+taille_bloc-1);
            blocs(:, idx) = b(:);
            idx = idx + 1;
        end
    end
end

function I = recombine_blocs(blocs, taille_bloc, L, C)
    I = zeros(L, C);
    idx = 1;
    for l = 1 : taille_bloc : L
        for c = 1 : taille_bloc : C
            I(l:l+taille_bloc-1, c:c+taille_bloc-1) = reshape(blocs(:, idx), taille_bloc, taille_bloc);
            idx = idx + 1;
        end
    end
end

% LBG avec splitting (section 10.4.1)
function [dictionnaire, D] = LBG_splitting(training_set, M)
    epsilon = 0.01;   % perturbation pour le splitting
    seuil   = 1e-3;   % seuil de convergence

    dictionnaire = mean(training_set, 2);
    while size(dictionnaire, 2) < M
        dictionnaire = [dictionnaire.*(1+epsilon), dictionnaire.*(1-epsilon)];
        D_prec = inf;
        while true
            [index, D] = encode_VQ(training_set, dictionnaire);
            for i = 1 : size(dictionnaire, 2)
                membres = training_set(:, index == i);
                if ~isempty(membres)          % "empty cell" (10.4.2) : on garde le codeword
                    dictionnaire(:, i) = mean(membres, 2);
                end
            end
            if D == 0 || abs(D_prec - D)/D < seuil
                break
            end
            D_prec = D;
        end
    end
end

function [index, D] = encode_VQ(vecteurs, dictionnaire)
    n = size(vecteurs, 2);
    index = zeros(1, n);
    dist_totale = 0;
    for k = 1 : n
        distances = sum((dictionnaire - vecteurs(:, k)).^2, 1);
        [d_min, index(k)] = min(distances);
        dist_totale = dist_totale + d_min;
    end
    D = dist_totale / n;
end

%% --- Quantificateur (sans variables globales) --------------------------
function Q = genere_quantif(nbits, Delta)
    N = 2^nbits;
    Q.codes = 0:(N-1);
    Q.seuils_decision = [-inf, ((-N/2+1):(N/2-1))*Delta, inf];
    Q.niveau_reconstruction = ((-N+1):2:(N-1)) * Delta/2;
end

function Out = Q93(Entree, Q)
    Out = zeros(size(Entree));
    for k = 1 : numel(Entree)
        x = Entree(k);
        for i = 1 : length(Q.codes)
            if x > Q.seuils_decision(i) && x <= Q.seuils_decision(i+1)
                Out(k) = Q.codes(i);
                break
            end
        end
    end
end

function Out = Q93_1(Entree, Q)
    Out = zeros(size(Entree));
    for k = 1 : numel(Entree)
        i = find(Q.codes == Entree(k), 1);
        Out(k) = Q.niveau_reconstruction(i);
    end
end

%% --- Utilitaires -------------------------------------------------------
function PSNR = calcul_PSNR(Iref, Idec, vmax)
    if nargin < 3
        vmax = 255;
    end
    eqm  = mean((double(Iref(:)) - double(Idec(:))).^2);
    PSNR = 10*log10(vmax^2 / eqm);
end

function afficher_image(num, I, titre)
    figure(num)
    imshow(uint8(I))
    title(titre)
end

function s = taille_str(I)
    s = size(I,1) + "x" + size(I,2);
end