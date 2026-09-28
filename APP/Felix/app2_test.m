%% Clear and load
clc; clear; close all;

Isource = imread("DPCM.bmp");

%% Problematique (Execution)
tic;

clc; close all;

% 1. Initialisation des parametres
L = 256;
C = 256;
N = 32; % Nombre de niveaux de quantification

% 2. Preparation de l'image
Iformat = formattageBiLin(Isource, L, C);
[Inorm, m, sigma] = normal(Iformat);

% 3. Recuperation des deltas theoriques (pour comparaison)
delta_base_quant = getDeltaUniforme(N);
delta_base_dpcm  = getDeltaLaplacien(N);

% 4. Plage globale d'optimisation de delta (balayage complet)
deltas_vecteur = linspace(1.5, 0.1, 100);

% 5. Boucles d'optimisation sur toute la plage de delta
psnr_quant_vals = zeros(size(deltas_vecteur));
psnr_dpcm_vals  = zeros(size(deltas_vecteur));

for i = 1:length(deltas_vecteur)
    d = deltas_vecteur(i);
    
    % Quantification Uniforme
    Iencode = Q93encoder_customDelta(Inorm, N, d);
    Idecode = Q93decoder_customDelta(Iencode, N, d);
    Idenorm = denormal(Idecode, m, sigma);
    psnr_quant_vals(i) = PSNR(Iformat, Idenorm);
    
    % DPCM
    Idpcm = DPCM2_customDelta(Iformat, N, d);
    psnr_dpcm_vals(i) = PSNR(Iformat, Idpcm);
end

% 6. Extraction des optimums globaux et valeurs de base
[max_psnr_quant, idx_opt_quant] = max(psnr_quant_vals);
delta_opt_quant = deltas_vecteur(idx_opt_quant);

[max_psnr_dpcm, idx_opt_dpcm] = max(psnr_dpcm_vals);
delta_opt_dpcm = deltas_vecteur(idx_opt_dpcm);

% PSNR avec les deltas theoriques de base
Iencode_base = Q93encoder_customDelta(Inorm, N, delta_base_quant);
Idecode_base = Q93decoder_customDelta(Iencode_base, N, delta_base_quant);
psnr_base_quant = PSNR(Iformat, denormal(Idecode_base, m, sigma));

psnr_base_dpcm = PSNR(Iformat, DPCM2_customDelta(Iformat, N, delta_base_dpcm));

% Reconstruction des images optimales
Idpcm_opt = DPCM2_customDelta(Iformat, N, delta_opt_dpcm);
Iencode_opt = Q93encoder_customDelta(Inorm, N, delta_opt_quant);
Idecode_opt = Q93decoder_customDelta(Iencode_opt, N, delta_opt_quant);
Idenorm_opt = denormal(Idecode_opt, m, sigma);

% 7. Affichage des resultats
figure('Name', sprintf('Optimisation Globale de Delta (N=%d)', N));

subplot(2,2,1);
plot(deltas_vecteur, psnr_quant_vals, 'b-', 'LineWidth', 1.5); hold on;
plot(delta_base_quant, psnr_base_quant, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k', ...
    'DisplayName', sprintf('\\delta_{théorie} = %.3f (%.2f dB)', delta_base_quant, psnr_base_quant));
plot(delta_opt_quant, max_psnr_quant, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', ...
    'DisplayName', sprintf('\\delta_{opt} = %.3f (%.2f dB)', delta_opt_quant, max_psnr_quant));
grid on; xlabel('\Delta'); ylabel('PSNR (dB)'); title('Quantification Uniforme'); legend('Location', 'best');

subplot(2,2,2);
plot(deltas_vecteur, psnr_dpcm_vals, 'r-', 'LineWidth', 1.5); hold on;
plot(delta_base_dpcm, psnr_base_dpcm, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k', ...
    'DisplayName', sprintf('\\delta_{théorie} = %.3f (%.2f dB)', delta_base_dpcm, psnr_base_dpcm));
plot(delta_opt_dpcm, max_psnr_dpcm, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
    'DisplayName', sprintf('\\delta_{opt} = %.3f (%.2f dB)', delta_opt_dpcm, max_psnr_dpcm));
grid on; xlabel('\Delta'); ylabel('PSNR (dB)'); title('DPCM (Loi Laplacienne)'); legend('Location', 'best');

subplot(2,2,3);
imshow(uint8(Idenorm_opt));
title(sprintf('Quant Opt (\\delta=%.3f) - PSNR: %.2f dB', delta_opt_quant, max_psnr_quant));

subplot(2,2,4);
imshow(uint8(Idpcm_opt));
title(sprintf('DPCM Opt (\\delta=%.3f) - PSNR: %.2f dB', delta_opt_dpcm, max_psnr_dpcm));

toc;


%% Fonctions du Traitement DPCM et Quantification

function [Idpcm, e, e_hat] = DPCM2_customDelta(Image, N, delta)
    Iint = convertToGray(Image);
    [L, C] = size(Iint);
    
    Idpcm = zeros(L, C);
    e = zeros(L, C);
    e_hat = zeros(L, C);
 
    sigma_e = std(Iint(:)); 
    step_size = delta * sigma_e;
    
    max_idx = (N / 2) - 1;
    min_idx = -(N / 2);

    for l = 1:L
        for c = 1:C
            p = pix_dpcm(Idpcm, l, c);
            
            e(l,c) = Iint(l,c) - p;
            e_quant_idx = min(max(round(e(l,c) / step_size), min_idx), max_idx);
            
            e_hat(l,c) = e_quant_idx * step_size;
            Idpcm(l,c) = min(max(p + e_hat(l,c), 0), 255);
        end
    end
end

function Iencode = Q93encoder_customDelta(Inorm, N, delta)
    seuils_int = ((-N/2 + 1):(N/2 - 1)) * delta;
    seuils = [-inf, seuils_int, inf];
    Iencode = discretize(Inorm, seuils) - 1;
end

function Idecode = Q93decoder_customDelta(Iencode, N, delta)
    seuils = ((-N + 1):2:(N - 1)) * delta / 2;
    Idecode = seuils(Iencode + 1);
end

function delta = getDeltaUniforme(N)
    nbits_dict = [2; 4; 6; 8; 10; 12; 14; 16; 32];
    delta_dict = [1.732; 0.866; 0.577; 0.433; 0.346; 0.289; 0.247; 0.217; 0.108];
    delta = delta_dict(nbits_dict == N);
end

function delta = getDeltaLaplacien(N)
    nbits_dict = [2; 4; 6; 8; 10; 12; 14; 16; 32];
    delta_dict = [1.414; 1.0873; 0.8707; 0.7309; 0.6334; 0.5613; 0.5055; 0.4609; 0.2799];
    delta = delta_dict(nbits_dict == N);
end

function psnr_val = PSNR(Isource, Imod)
    sigma2 = mean((double(Isource(:)) - double(Imod(:))).^2);
    psnr_val = 10 * log10(max(double(Isource(:)))^2 / sigma2);
end

function Iformat = formattageBiLin(Isource, L2, C2)
    Iint = convertToGray(Isource);
    [L1, C1] = size(Iint);
    
    r = (1:L2)' * (L1 / L2) - 0.5 * (L1 / L2 - 1);
    k = (1:C2) * (C1 / C2) - 0.5 * (C1 / C2 - 1);
    
    ll = min(max(floor(r), 1), L1 - 1);
    cc = min(max(floor(k), 1), C1 - 1);
    
    dl = r - ll;
    dc = k - cc;
    
    A = Iint(ll, cc);
    B = Iint(ll, cc + 1);
    C = Iint(ll + 1, cc);
    D = Iint(ll + 1, cc + 1);
    
    X = (1 - dc) .* A + dc .* B;
    Y = (1 - dc) .* C + dc .* D;
    
    Iformat = floor((1 - dl) .* X + dl .* Y);
end

function [Inorm, m, sigma] = normal(Image)
    doubleImg = double(Image);
    m = mean(doubleImg(:));
    sigma = std(doubleImg(:));
    Inorm = (doubleImg - m) / sigma;
end

function Idenorm = denormal(Image, m, sigma)
    Idenorm = (Image * sigma) + m;
end

function p = pix_dpcm(p_pred, L, C)
    if L == 1 && C == 1
        p = 128;
    elseif L == 1
        p = p_pred(1, C - 1);
    elseif C == 1
        p = p_pred(L - 1, 1);
    else
        p1 = 0.5 * p_pred(L-1, C)   + 0.5 * p_pred(L, C-1);
        p2 = 0.5 * p_pred(L-1, C-1) + 0.5 * p_pred(L, C-1);
        p3 = 0.5 * p_pred(L-1, C-1) + 0.5 * p_pred(L-1, C);
        p = median([p1, p2, p3]);
    end
end

function Igray = convertToGray(Image)
    if size(Image, 3) == 3
        Igray = mean(double(Image), 3);
    else
        Igray = double(Image);
    end
end