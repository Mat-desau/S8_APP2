%% Load and clear



%% Q1

clc; clear all; close all;
image = imread("cman.tif");
[L1, C1] = size(image);

L2 = 128;
C2 = 128;

Iredim = zeros(L2, C2);

for l=1 : L2
    for c=1 : C2
        ll = floor(l*L1/L2);
        cc = floor(c*C1/C2);
        Iredim(l,c) = image(ll,cc);
    end
end

imshow(uint8(Iredim))

%% Q2

clc; clear all; close all;
image = imread("lenna.bmp");
[L1, C1, Z1] = size(image);

R = double(image(:,:,1));
G = double(image(:,:,2));
B = double(image(:,:,3));

imod = (R+G+B)/3;

L2 = 256;
C2 = 256;


Iredim = zeros(L2, C2);

for l=1 : L2
    for c=1 : C2
        ll = floor(l*L1/L2);
        cc = floor(c*C1/C2);
        Iredim(l,c) = imod(ll,cc);
    end
end

imshow(uint8(Iredim))

%% Q3

clc; clear all; close all;
Isource = imread("cman.tif");
[L1, C1] = size(Isource);

m = mean(Isource(:));
sigma = std(double(Isource(:)));

delta = 0.433;
n = 8;

Inorm = (double(Isource)-m)/sigma;

Iencode = Q93encoder(Inorm, delta, n);

Idecode = Q93decoder(Iencode, delta, n);

Ifin = (Idecode*sigma)+m;

imshow(Ifin);

%% Fonctions

function encode = Q93encoder(Inorm, delta, n)
    seuils = [-inf, -3*delta, -2*delta, -delta, 0 delta, 2*delta, 3*delta, inf];
    codes = [0, 1, 2, 3, 4, 5, 6, 7];

    [L,C] = size(Inorm);
    encode = zeros(L,C);

    for l = 1:L
        for c = 1:C
            x = Inorm(l,c);
            for i = 1:n
                if x > seuils(i) && x <= seuils(i+1)
                    encode(l,c) = codes(i);
                    break;
                end
            end
        end
    end
end

function decode = Q93decoder(Iencode, delta, n)
    seuils = [-7*delta/2, -5*delta/2, -3*delta/2, -delta/2, delta/2, 3*delta/2, 5*delta/2, 7*delta/2];
    codes = [0, 1, 2, 3, 4, 5, 6, 7];

    [L,C] = size(Iencode);
    decode = zeros(L,C);

    for l = 1:L
        for c = 1:C 
            for i = 1:n
                if Iencode(l,c) == codes(i);
                    decode(l,c) = seuils(i);
                    break;
                end
            end
        end
    end
end