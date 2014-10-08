clear all; close all;

load('DataHW3.mat')

% Remove outliers

ind = P < 60;
P = P(ind);
T = T(ind);
E = E(ind);
N = N(ind);

%% Problem 1
% (3 points)
% (a)

meanT = mean(T);
disp(['mean of T: ', num2str(meanT)])
meanP = mean(P);
disp(['mean of P: ', num2str(meanP)])

stdT = std(T);
disp(['std of T: ', num2str(stdT)])
stdP = std(P);
disp(['std of P: ', num2str(stdP)])

%%
% (b)

step = 2;
bounds = 0:step:10;
nw = floor(10/step);
Tw = nan(nw);
Pw = nan(nw);


for i=1:nw
    for j=1:nw
        ind = E>bounds(j) & E <= bounds(j+1) & ...
              N>bounds(i) & N <= bounds(i+1);
        Tw(i, j) = mean(T(ind));
        Pw(i, j) = mean(P(ind));
    end
end

meanT_glob = mean2(Tw);
meanP_glob = mean2(Pw);

disp(['Global mean of T (2x2): ', num2str(meanT_glob)])
disp(['Global mean of P (2x2): ', num2str(meanP_glob)])

%%
% (c)

Wsize = .5:.5:10;
n = length(Wsize);
meanT = zeros(n, 1);
meanP = zeros(n, 1);


for k=1:n
    step = Wsize(k);
    bounds = 0:step:10;
    nw = floor(10/step);
    Tw = nan(nw);
    Pw = nan(nw);


    for i=1:nw
        for j=1:nw
            ind = E>bounds(j) & E <= bounds(j+1) & ...
                  N>bounds(i) & N <= bounds(i+1);
            Tw(i, j) = mean(T(ind));
            Pw(i, j) = mean(P(ind));
        end
    end
    
    meanT(k) = mean2(Tw);
    meanP(k) = mean2(Pw);
end

figure
plot(Wsize, meanT)

figure
plot(Wsize, meanP)

%% Problem 2
% (4 points)
% (a)
E1 = 1; N1 = 6;
E2 = 6; N2 = 4;

%%
% polygonal method
n = length(N);
dist1 = Inf; dist2 = Inf;
P1 = NaN; P2 = NaN; T1 = NaN; T2 = NaN; 
for i=1:n
    dN1 = N1 - N(i);
    dE1 = E1 - E(i);
    dist1_tmp = sqrt(dN1^2 + dE1^2);
    if (dist1_tmp < dist1)
        dist1 = dist1_tmp;
        P1 = P(i); T1 = T(i);
    end
    
    dN2 = N2 - N(i);
    dE2 = E2 - E(i);
    dist2_tmp = sqrt(dN2^2 + dE2^2);
    if (dist2_tmp < dist2)
        dist2 = dist2_tmp;
        P2 = P(i); T2 = T(i);
    end
end

disp(['P(x1, y1) = ' , num2str(P1)])
disp(['P(x2, y2) = ' , num2str(P2)])

disp(['T(x1, y1) = ' , num2str(T1)])
disp(['T(x2, y2) = ' , num2str(T2)])

%%
% triangular method
Pinterp = TriScatteredInterp(E, N, P);
P1 = Pinterp(E1, N1);
P2 = Pinterp(E2, N2);

Tinterp = TriScatteredInterp(E, N, T);
T1 = Tinterp(E1, N1);
T2 = Tinterp(E2, N2);

disp(['P(x1, y1) = ' , num2str(P1)])
disp(['P(x2, y2) = ' , num2str(P2)])

disp(['T(x1, y1) = ' , num2str(T1)])
disp(['T(x2, y2) = ' , num2str(T2)])

%%
% Inverse distance method

dN1 = N - N1;
DE1 = E - E1;
dist1 = sqrt(dN1.^2 + dE1.^2);
P1 = dist1' * P / sum(dist1); 
T1 = dist1' * T / sum(dist1);

dN2 = N - N2;
dE2 = E - E2;
dist2 = sqrt(dN2.^2 + dE2.^2);
P2 = dist2' * P / sum(dist2);
T2 = dist2' * T / sum(dist2);

disp(['P(x1, y1) = ' , num2str(P1)])
disp(['P(x2, y2) = ' , num2str(P2)])

disp(['T(x1, y1) = ' , num2str(T1)])
disp(['T(x2, y2) = ' , num2str(T2)])


%%
% (b)

spacing = 0:.1:10;
n = length(spacing);
Tgrid = nan(n);

for i=1:n
    for j = 1:n
        Ns = spacing(i);
        Es = spacing(j);
        dN = Ns - N;
        dE = Es - E;
        dist = sqrt(dN.^2 + dE.^2);
        Tgrid(i, j) = dist' * T / sum(dist);
    end
end

imagesc(Tgrid)

%%
spacing = 0:.1:10;
n = length(spacing);
Tgrid = nan(n);

for i=1:n
    for j = 1:n
        Ns = spacing(i);
        Es = spacing(j);
        Tgrid(i, j) = Tinterp(Es, Ns);
    end
end

imagesc(Tgrid)
