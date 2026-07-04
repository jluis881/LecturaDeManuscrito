clc;
clear;
close all;

%% ===============================
% NUMEROS MANUSCRITOS
% EULER + CHAIN + CORR
% ===============================


%% ===============================
% CARGA DE PLANTILLAS
% ===============================
    
    templates = cell(10,1);
    chainTemplates = cell(10,1);
    
for n=0:9
    img = imread(sprintf('%d.jpg',n));
    img = rgb2gray(img);
    img = adapthisteq(img);
    img = imbinarize(img);
    img = ~img;
    img = bwareaopen(img,20);
    img = imclose(img,strel('disk',2));
    img = preprocess(img);
  
    templates{n+1} = img;
    chainTemplates{n+1} = slopeDescriptor(img);
end

%% ===============================
% ANALISIS DE LAS IMAGENES
% ===============================

for img=1:4
    I = imread(sprintf('numeros%d.jpg',img));
    I = rgb2gray(I);
    I = adapthisteq(I);
    
    BW = imbinarize(I);
    BW = ~BW;
    BW = bwareaopen(BW,40);
    BW = imclose(BW,strel('disk',2)); %Usar dilatación si hay números abiertos
    
    [L,~] = bwlabel(BW,8);
    
    props = regionprops(L,'Image','BoundingBox','Area');
    props = props([props.Area] > 80);
    
    x = arrayfun(@(p)p.BoundingBox(1),props);
    [~,idx] = sort(x);
    props = props(idx);
    
    %% ===============================
    % MOSTRAR DETECCION
    % ===============================
    
    figure;
    imshow(I);
    hold on;
    
    for k=1:length(props)
        rectangle('Position',props(k).BoundingBox,...
            'EdgeColor','r','LineWidth',2);
    end
    title(sprintf('Digitos detectados de imagen %d',img));
     
    %% ===============================
    % CLASIFICACION
    % ===============================
    
    codigo = '';
    figure;
    for k=1:length(props)
        digit = preprocess(props(k).Image);
    
        %% =========================
        % ETAPA 1: EULER
        %% =========================
    
        E = bweuler(digit);
        if E == -1
            candidatos = 8;
        elseif E == 0
            candidatos = [0 6 9];
        else
            candidatos = [1 2 3 4 5 7];
        end
    
        %% =========================
        % ETAPA 2: CHAIN 
        %% =========================
    
        desc = slopeDescriptor(digit);
        bestChainDist = inf;
        reducedCandidates = [];
    
        for c = candidatos
            d = norm(desc - chainTemplates{c+1});
            if d < bestChainDist + 0.05   % tolerancia ligera
                reducedCandidates(end+1) = c;
            end
        end
        if isempty(reducedCandidates)
            reducedCandidates = candidatos;
        end
    
        %% =========================
        % ETAPA 3: CORR2 FINAL
        %% =========================
    
        bestScore = -inf;
        bestDigit = '?';
        for c = reducedCandidates
            score = corr2_sliding(digit, templates{c+1});
            if score > bestScore
                bestScore = score;
                bestDigit = num2str(c);
            end
        end
    
        codigo = [codigo bestDigit];
        subplot(ceil(length(props)/5),5,k);
        imshow(digit);
        title(bestDigit);
    end
    
    disp('=====================');
    disp('CODIGO DETECTADO');
    disp(codigo);
    
    % Evaluar coincidencia de las imagenes
    objetivo = '7890345645123456789';
    if ~ischar(codigo) && ~isstring(codigo)
        codigo = char(codigo);
    else
        codigo = char(codigo);
    end
    nComp = min(length(codigo), length(objetivo));
    if nComp == 0
        porcentaje = 0;
    else
        iguales = sum(codigo(1:nComp) == objetivo(1:nComp));
        porcentaje = iguales / length(objetivo) * 100; % porcentaje respecto a toda la cadena objetivo
    end
    % Mostrar resultados
    fprintf('Codigo detectado: %s\n', codigo);
    fprintf('Cadena objetivo:  %s\n', objetivo);
    fprintf('Coincidencias en los primeros %d caracteres: %d\n', nComp, iguales);
    fprintf('Porcentaje de coincidencia: %.2f%%\n', porcentaje);
    disp('=====================');
end



%% ===============================
% PREPROCESADO
% ===============================

function out = preprocess(img)

    img = logical(img);
    [r,c] = find(img);
    if isempty(r)
        out = zeros(60,40);
        return;
    end

    img = img(min(r):max(r),min(c):max(c));
    img = imresize(img,[50 30]);
    out = zeros(60,40);
    r0 = floor((60-size(img,1))/2);
    c0 = floor((40-size(img,2))/2);
    out(r0+1:r0+size(img,1),c0+1:c0+size(img,2)) = img;
end

%% ===============================
% CHAIN DESCRIPTOR
% ===============================

function desc = slopeDescriptor(img)

    B = bwboundaries(img,'noholes');
    if isempty(B)
        desc = zeros(1,8);
        return;
    end

    P = B{1};
    dirs = zeros(size(P,1)-1,1);

    for k=1:size(P,1)-1
        dy = P(k+1,1)-P(k,1);
        dx = P(k+1,2)-P(k,2);
        ang = atan2d(dy,dx);
        if ang < 0
            ang = ang + 360;
        end
        dirs(k) = floor(ang/45);
    end

    diffDirs = mod(diff(dirs),8);
    desc = histcounts(diffDirs,-0.5:7.5);
    desc = desc ./ (sum(desc)+eps);
end

%% ===============================
% FUNCION CORRELACION
% ===============================

function score = corr2_sliding(A,B)

    A = double(A);
    B = double(B);
    best = -inf;

    for i=-3:3
        for j=-3:3
            Bshift = circshift(B,[i j]);
            c = corr2(A,Bshift);

            if c > best
                best = c;
            end
        end
    end

    score = best;
end
