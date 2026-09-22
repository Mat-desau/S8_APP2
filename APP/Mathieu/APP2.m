%% APP2
% BOIF1302
% DESM1210

clc
clear all
close all

%% Variable global
global seuils_decision codes niveau_reconstruction

%% Problematique
% Loader les images
% Lenna = 512x512
% cman = 256x256
% irm = 256x256
% mandrill = 256x256


% Isource = imread("lenna.bmp");
% Isource = imread("cman.tif");
Isource = imread("irm.tif");
% Isource = imread("mandrill.tif");
% Isource = imread("crest.bmp");
[L1, C1, Z] = size(Isource);

if Z == 3
    % Pour avoir le RGB 
    R = double(Isource(:,:,1));
    G = double(Isource(:,:,2));
    B = double(Isource(:,:,3));
    % Refaire la modifier en 1 couleur
    IsourceMod = (R+G+B)/3;
else
    IsourceMod = double(Isource);
end

%% Redimensionnement PPV
% Analyse si c'est pas déjà 256
if L1 > 256 && C1 > 256
    L2 = 256;
    C2 = 256;
    Iredim = zeros(L2, C2);
    
    % FOR pour redimensioner
    for l = 1 : L2
        for c = 1 : C2
            ll = floor(l*L1/L2);
            cc = floor(c*C1/C2);
            if ll == 0
                ll = 1;
            end
            if cc == 0
                cc = 1;
            end
            Iredim_PPV(l, c) = IsourceMod(ll, cc);
        end
    end
else
    Iredim_PPV = IsourceMod;
    L2 = L1;
    C2 = C1;
end

%% Redimensionnement Bilineaire
if L1 > 256 && C1 > 256
    L2 = 256;
    C2 = 256;
    Iredim = zeros(L2, C2);
    
    % FOR pour redimensioner
    for l = 1 : L2-1
        for c = 1 : C2-1
            ll = floor(l*L1/L2);
            cc = floor(c*C1/C2);
            if ll == 0
                ll = 1;
            end
            if cc == 0
                cc = 1;
            end
            A = IsourceMod(ll, cc);
            B = IsourceMod(ll, cc+1);
            C = IsourceMod(ll+1, cc);
            D = IsourceMod(ll+1, cc+1);
            X = (1-(c*C1/C2-cc))*A+(c*C1/C2-cc)*B;
            Y = (1-(c*C1/C2-cc))*C+(c*C1/C2-cc)*D;
            Iredim_BiL(l, c) = uint8(floor((1-(l*L1/L2-ll))*X+(l*L1/L2-ll)*Y));
        end
    end
else
    Iredim_BiL = IsourceMod;
    L2 = L1;
    C2 = C1;
end

%% Quantificateur Scalaire (QS)
m = mean(Iredim_PPV(:));
sigma = std(double(Iredim_PPV(:)));

% Bits = [2, 4, 6, 8, 10, 12, 14, 16, 32];
% Delta_vals = [1.414, 1.0873, 0.8707, 0.7309, 0.6334, 0.5613, 0.5055, 0.4609, 0.2799];
% Delta_vals = [1.596, 0.9957, 0.7334, 0.5860, 0.4908, 0.4238, 0.3739, 0.3352, 0.1881];
% delta_map = containers.Map(Bits, Delta_vals);

nbits = 5;
N = 2^nbits;
Delta = 0.1881;

[seuils_decision, codes, niveau_reconstruction] = genere_quantif(nbits, Delta);

for l = 1 : L2
    for c = 1 : C2
        Inormaliser(l, c) = (Iredim_PPV(l, c)-m)/sigma;
        Icode(l, c) = Q93(Inormaliser(l,c));
        InormQ(l, c) = Q93_1(Icode(l, c));
        IDecoder_QS(l, c) = (InormQ(l, c)*sigma)+m;
    end
end
sigma2 = mean((double(Iredim_PPV(:))-double(IDecoder_QS(:))).^2);
PSNR_QS = 10 * log10((255^2)/(sigma2));

%% Quantificateur Différentiel (DPCM)
Icode_DPCM = zeros(L2, C2);

nbits = 5;
N = 2^nbits;
Delta = 0.2779;

[seuils_decision, codes, niveau_reconstruction] = genere_quantif(nbits, Delta);


for l = 1 : L2
    for c = 1 : C2
        % Prédiction
        if l == 1 && c == 1
            xpred = 128;
        elseif l == 1
            xpred = Iredim_PPV(l, c-1);
        elseif c == 1
            xpred = Iredim_PPV(l-1, c);
        else
            p1 = 0.5*Iredim_PPV(l-1, c)   + 0.5*Iredim_PPV(l, c-1);
            p2 = 0.5*Iredim_PPV(l-1, c-1) + 0.5*Iredim_PPV(l, c-1);
            p3 = 0.5*Iredim_PPV(l-1, c-1) + 0.5*Iredim_PPV(l-1, c);
            xpred = median([p1, p2, p3]);
        end

        % Erreur de prédiction
        erreur_DPCM(l,c) = Iredim_PPV(l, c) - xpred;
    end
end

m = mean(erreur_DPCM(:));
sigma = std(double(erreur_DPCM(:)));

for l = 1 : L2
    for c = 1 : C2
        % Prédiction
        if l == 1 && c == 1
            xpred = 128;
        elseif l == 1
            xpred = Iredim_PPV(l, c-1);
        elseif c == 1
            xpred = Iredim_PPV(l-1, c);
        else
            p1 = 0.5*Iredim_PPV(l-1, c)   + 0.5*Iredim_PPV(l, c-1);
            p2 = 0.5*Iredim_PPV(l-1, c-1) + 0.5*Iredim_PPV(l, c-1);
            p3 = 0.5*Iredim_PPV(l-1, c-1) + 0.5*Iredim_PPV(l-1, c);
            xpred = median([p1, p2, p3]);
        end

        % Erreur de prédiction
        e = Iredim_PPV(l, c) - xpred;

        I_DPCM = (e-m)/sigma;

        % Quantification de l'erreur
        Icode_DPCM = Q93(I_DPCM);
        eq = Q93_1(Icode_DPCM);
    
        % Reconstruction
        IDecoder_DPCM(l, c) = ((eq*sigma)+m)+xpred;
    end
end

sigma2 = mean((double(Iredim_PPV(:))-double(IDecoder_DPCM(:))).^2);
PSNR_DPCM = 10 * log10((255^2)/(sigma2));

%% Affichage
figure(1)
imshow(Isource);
title("ISource " + L1 + "x" + C1)
figure(2)
imshow(uint8(IsourceMod))
title("ISourceMod " + L1 + "x" + C1)
figure(3)
imshow(uint8(Iredim_PPV))
title("IRedim PPV" + L2 + "x" + C2)
figure(4)
imshow(uint8(Iredim_BiL))
title("IRedim BiL" + L2 + "x" + C2)
figure(5)
imshow(uint8(IDecoder_QS))
title("IDecoder QS " + N + " niveau PSNR = " + PSNR_QS)
figure(6)
imshow(uint8(IDecoder_DPCM))
title("IDecoder DPCM " + N + " niveau PSNR = " + PSNR_DPCM)


%% Fonction Q93
function Icode = Q93(Inormaliser)
    global seuils_decision codes

    [L, C] = size(Inormaliser);
    Icode = zeros(L, C);

    for l = 1:L
        for c = 1:C
            x = Inormaliser(l, c);
            for i = 1:length(codes)
                if x > seuils_decision(i) && x <= seuils_decision(i+1)
                    Icode(l, c) = codes(i);
                    break;
                end
            end
        end
    end
end

%% Fonction Q93_1
function InormQ = Q93_1(Icode)
    global codes niveau_reconstruction

    [L, C] = size(Icode);
    InormQ = zeros(L, C);

    for l = 1:L
        for c = 1:C
            for i = 1 : length(codes)
                if Icode(l,c) == codes(i)
                    InormQ(l, c) = niveau_reconstruction(i);
                    break;
                end
            end
        end
    end
end

%% Generation de quantification
function [seuils_decision, codes, niveau_reconstruction] = genere_quantif(nbits, Delta)
    N = 2^nbits;  % nombre de niveaux, ex: 8 pour nbits=3
    codes = 0:(N-1);

    seuils_internes = ((-N/2+1):(N/2-1)) * Delta;
    seuils_decision = [-inf, seuils_internes, inf];

    niveau_reconstruction = ((-N+1):2:(N-1)) * Delta/2;
end