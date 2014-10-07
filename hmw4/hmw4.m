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


%%
% (b)
