close all;
clear all;

%% Problem 1
% (3 points)

load DataHW3.mat

%%
% a) (0.5 points)

disp(['Mean of T: ' num2str(mean2(T))])
disp(['Mean of P: ' num2str(mean2(P))])

disp(['Median of T: ' num2str(median(T(:)))])
disp(['Median of P: ' num2str(median(P(:)))])

%% 
% b) (0.5 points)

figure
hist(T, 20);
title('Histogram of T')
xlabel('T')
ylabel('Count')

disp(['Skewness of T: ', num2str(skewness(T(:)))])

% The histogram of T and the skewness coefficient show that the T is
% positively skewed.

%%

figure
hist(P, 20)
title('histogram of P with outliers')
ylabel('count')
xlabel('P')

% We notice that most of the values of P are less than 100. Because of a
% few very large outliers, we cannot really see  the distribution of most 
% of the values of P. 

%%

figure
hist(P(P < 60), 20);
title('P without outliers')
ylabel('count')
xlabel('P')

% Now we see that the distribution of P is monomodal, almost symmetric,
% with perhaps a slight negative skewness.

%%

disp(['Skewness of P with outlier: ', num2str(skewness(P(:)))])

% If we compute the skewness of P with all values we find a quite large
% positive value. However, we know that the skewness coefficient is very
% sensitive to outliers. 

%% 

disp(['Skewness of P without outliers: ', ...
    num2str(skewness(P(P < 60)))])

% If we remove a few extreme values we find a small negative values which
% agrees with what we observe in the histogram.

%%
% c) (1 point)


% To see if the P and T are correlated we first plot a scatter plot:

scatter(T, P)
xlabel('T')
ylabel('P')
title('Scatter plot of P and T with outliers')

% Again, because of some extreme values taken by P, we can't really see how
% most of the values behave. 

%%

scatter(P(P<100), T(P<100))
xlabel('T')
ylabel('P')
title('Scatter plot of P and T without outliers')

% Now we can see a negative correlation. 

%%
% We try to confirm the negative correlation by computing the correlation
% coefficient.

disp(['Correlation coeff with outliers: ', num2str(corr(P, T))])

% This is a weak correlation. But it may be influenced by the outliers.
% Indeed if we plot the regression line with this coefficient we can see
% that it does not really fit the data:


disp(['Correlation coeff without outliers: ',...
    num2str(corr(P(P<100), T(P<100)))])

% This corresponds better to what we observe. 

p1 = polyfit(P, T, 1);
p2 = polyfit(P(P<100), T(P<100), 1);
figure
scatter(P(P<100), T(P<100))
hold on
plot(0:80, polyval(p1, 0:80), 'r')
plot(0:80, polyval(p2, 0:80), 'k')
xlabel('P')
ylabel('T')
title('scatter plot with regression lines')
legend('observations', 'with outliers', 'without outliers')

%%
% d)

% From the histograms in b), we clearly see that P has some very large
% isolated values, that do not seem to be distributed like the rest of the 
% values. T has some very large values too but they fit the general
% skewness of the distribution. In c) we showed that the large majority of
% the data values of P tend to decrease for increasing values of T. The
% 4 value of P above 60 are both isolated from the other values and do not
% fit this trend. We decide to flag them as outliers. T has some large
% values too, but they fit the decreasing trend, so we keep them.

figure
scatter(P(P<60), T(P<60), 'b')
hold on
scatter(P(P>60), T(P>60), 'r')
plot(0:500, polyval(p1, 0:500), 'r')
xlabel('P')
ylabel('T')
title('scatter plot with outliers')
legend('"good" values', 'outliers')

%% Problem 2
% (7 points)

load DataHW3.mat

%%
% a) (0.5 points)

plot(E, N, 'xk')
xlabel('Easting')
ylabel('Northing')
title('Position of the observations')

% We can see that the data the sparsest distributions in the northern and
% north-western areas. Its is most densely distributed along curves that
% are distributed over the area.

%%
% b) (0.5 points)

figure
scatter(E, N, 50, T, 'filled')
colorbar
xlabel('Easting')
ylabel('Northing')
title('Values of T')


figure
scatter(E, N, 50, P, 'filled')
caxis([0, 60])
colorbar
xlabel('Easting')
ylabel('Northing')
title('Values of P')

% caxis function will force to scale the colors so we don't only see the
% outliers


%% 
% c) (1 point)
% We divide the area into 10*10 cells. We can clearly split the easting and
% the northing by ten between 0 and 10. 

n = 10;
% 
% Emin = min(E); Emax = max(E); Erange = Emax - Emin; dE = Erange/n;
% Nmin = min(N); Nmax = max(N); Nrange = Nmax - Nmin; dN = Nrange/n;
% 
% Ebounds = Emin:dE:Emax;
% Nbounds = Nmin:dN:Nmax;

Ebounds = 0:10;
Nbounds = 0:10;

P_meanW = zeros(n); P_stdW  = zeros(n);

for i=1:n
    for j=1:n
        P_meanW(j,i) = mean(P(E >= Ebounds(i) & E < Ebounds(i+1) & ...
                              N >= Nbounds(j) & N < Nbounds(j+1)));
        P_stdW(j,i)  = std(P(E >= Ebounds(i) & E < Ebounds(i+1) & ...
                              N >= Nbounds(j) & N < Nbounds(j+1)));
    end
end

figure
imagesc(.5:1:9.5, .5:1:9.5, P_meanW)
xlabel('Easting')
ylabel('Northing')
colorbar
caxis([0, 60])
title('Moving-window average')
set(gca,'YDir','normal')

figure
imagesc(.5:1:9.5, .5:1:9.5, P_stdW)
xlabel('Easting')
ylabel('Northing')
colorbar
caxis([0, 20])
title('Moving-window std')
set(gca,'YDir','normal')

%%
% d) (1 point)


figure
scatter(P_meanW(:), P_stdW(:))
xlabel('means')
ylabel('std')
title('with outliers')
% The large outliers prevent us from seeing if there is trend for most of
% the data.

figure
scatter(P_meanW(:), P_stdW(:))
xlim([15 50])
ylim([0  15])
xlabel('means')
ylabel('std')
title('without outliers')
% When we zoom where most of the data points are, we see no major
% correlation.

% If we compute the correlation coefficient between the mean and the
% standard deviation we find:

corr(P_meanW(:), P_stdW(:))

% It is a quite important correlation but it might be caused by the
% outliers. When we exclude them, we find:

corr(P_meanW(P_stdW < 15), P_stdW(P_stdW < 15))

% This is a very low correlation. Therefore, we can't conclude that there 
% is a proportional effect in P.

%%
% e) (1 point)
% We bin distances from 0 to 10, with intervals of size 0.5. 

bins = [-Inf, 0:.5:10, +Inf];
nbBins = size(bins, 2) - 1;

[X, Y] = meshgrid(E, N);
distances = sqrt((X - X').^2 + (Y - Y').^2);
variog = zeros(nbBins, 1);
centers = zeros(nbBins, 1);

for i=1:nbBins
    [row, col] = find(distances > bins(i) & distances <= bins(i+1));
    centers(i) = mean(distances(sub2ind(size(distances), row, col)));
    variog(i) = 0.5 * mean((T(row) - T(col)).^2);
end

figure
plot(centers, variog)
title('Variogram of T')

%%
% Since the variogram is sensible to outlier, we remove them from P.

bins = [-Inf, 0:.5:10, +Inf];
nbBins = size(bins, 2) - 1;
ind = find(P < 80);
Eb = E(ind); Nb = N(ind); Pb = P(ind);

[X, Y] = meshgrid(Eb, Nb);
distances = sqrt((X - X').^2 + (Y - Y').^2);
variog = zeros(nbBins, 1);
centers = zeros(nbBins, 1);

for i=1:nbBins
    [row, col] = find(distances > bins(i) & distances <= bins(i+1));
    centers(i) = mean(distances(sub2ind(size(distances), row, col)));
    variog(i) = 0.5 * mean((Pb(row) - Pb(col)).^2);
end

figure
plot(centers, variog)
title('Variogram of P')


%%
% f) (1 point)
bins = [0:.5:9, +Inf];
nbBins = size(bins, 2) - 1;

xVariog = zeros(nbBins+1, 1);
xCenters = zeros(nbBins+1, 1);

yVariog = zeros(nbBins+1, 1);
yCenters = zeros(nbBins+1, 1);

n = length(T);
angtol = 15;

for k=1:nbBins
    xcount = 0;
    ycount = 0;
    x_tmp = 0;
    y_tmp = 0;
    x_cen = 0;
    y_cen = 0;
    for i=1:n
        for j=1:n
            if (i~=j)
                dN = N(i) - N(j); dE = E(i) - E(j);
                dist = sqrt(dN^2 + dE^2);
                ang = atan(dN/dE) * pi / 180;
                if (dist > bins(k) && dist <= bins(k+1))
                    if (ang < angtol && ang > -angtol)
                        x_tmp = x_tmp + (T(i) - T(j))^2;
                        xcount = xcount + 1;
                        x_cen = x_cen + dist;
                    end
                    if (dE == 0 || ang < -90 + angtol && ang > 90 - angtol)
                        y_tmp = y_tmp + (T(i) - T(j))^2;
                        ycount = ycount + 1;
                        y_cen = y_cen + dist;
                    end
                end
            end
        end
    end
    xVariog(k+1) = 0.5 * x_tmp / xcount;
    yVariog(k+1) = 0.5 * y_tmp / ycount;
    xCenters(k+1) = x_cen / xcount;
    yCenters(k+1) = y_cen / ycount;
end

        
% 
% for i=1:nbBins
%     [row, col] = find(xDistances > bins(i) & xDistances <= bins(i+1));
%     xCenters(i) = mean(xDistances(sub2ind(size(xDistances), row, col)));
%     xVariog(i) = 0.5 * mean((T(row) - T(col)).^2);
% end
% 
% for i=1:nbBins
%     [row, col] = find(yDistances > bins(i) & yDistances <= bins(i+1));
%     yCenters(i) = mean(yDistances(sub2ind(size(yDistances), row, col)));
%     yVariog(i) = 0.5 * mean((T(row) - T(col)).^2);
% end

figure
plot(yCenters, yVariog, '-or')
hold on
plot(xCenters, xVariog, '-ob')
hold off
title('Directional variograms of T')
legend('N-S', 'E-W')


%%

% Exclude outliers
Pc = P(P < 60);
Ec = E(P < 60);
Nc = N(P < 60);

bins = [0:.5:9, +Inf];
nbBins = size(bins, 2) - 1;

xVariog = zeros(nbBins+1, 1);
xCenters = zeros(nbBins+1, 1);

yVariog = zeros(nbBins+1, 1);
yCenters = zeros(nbBins+1, 1);

n = length(Pc);
angtol = 15;

for k=1:nbBins
    xcount = 0;
    ycount = 0;
    x_tmp = 0;
    y_tmp = 0;
    x_cen = 0;
    y_cen = 0;
    for i=1:n
        for j=1:n
            if (i~=j)
                dN = Nc(i) - Nc(j); dE = Ec(i) - Ec(j);
                dist = sqrt(dN^2 + dE^2);
                ang = atan(dN/dE) * pi / 180;
                if (dist > bins(k) && dist <= bins(k+1))
                    if (ang < angtol && ang > -angtol)
                        x_tmp = x_tmp + (Pc(i) - Pc(j))^2;
                        xcount = xcount + 1;
                        x_cen = x_cen + dist;
                    end
                    if (dE == 0 || ang < -90 + angtol && ang > 90 - angtol)
                        y_tmp = y_tmp + (Pc(i) - Pc(j))^2;
                        ycount = ycount + 1;
                        y_cen = y_cen + dist;
                    end
                end
            end
        end
    end
    xVariog(k+1) = 0.5 * x_tmp / xcount;
    yVariog(k+1) = 0.5 * y_tmp / ycount;
    xCenters(k+1) = x_cen / xcount;
    yCenters(k+1) = y_cen / ycount;
end

        
% 
% for i=1:nbBins
%     [row, col] = find(xDistances > bins(i) & xDistances <= bins(i+1));
%     xCenters(i) = mean(xDistances(sub2ind(size(xDistances), row, col)));
%     xVariog(i) = 0.5 * mean((T(row) - T(col)).^2);
% end
% 
% for i=1:nbBins
%     [row, col] = find(yDistances > bins(i) & yDistances <= bins(i+1));
%     yCenters(i) = mean(yDistances(sub2ind(size(yDistances), row, col)));
%     yVariog(i) = 0.5 * mean((T(row) - T(col)).^2);
% end

figure
plot(yCenters, yVariog, '-or')
hold on
plot(xCenters, xVariog, '-ob')
hold off
title('Directional variograms of P (without outliers)')
legend('N-S', 'E-W')

%%






bins = [-Inf, 0:.5:10, +Inf];
nbBins = size(bins, 2) - 1;
ind = find(P < 80);
Eb = E(ind); Nb = N(ind); Pb = P(ind);

[X, Y] = meshgrid(Eb, Nb);

xDistances = abs(X - X');
xVariog = zeros(nbBins, 1);
xCenters = zeros(nbBins, 1);

yDistances = abs(Y - Y');
yVariog = zeros(nbBins, 1);
yCenters = zeros(nbBins, 1);

for i=1:nbBins
    [row, col] = find(xDistances > bins(i) & xDistances <= bins(i+1));
    xCenters(i) = mean(xDistances(sub2ind(size(xDistances), row, col)));
    xVariog(i) = 0.5 * mean((Pb(row) - Pb(col)).^2);
end

for i=1:nbBins
    [row, col] = find(yDistances > bins(i) & yDistances <= bins(i+1));
    yCenters(i) = mean(yDistances(sub2ind(size(yDistances), row, col)));
    yVariog(i) = 0.5 * mean((Pb(row) - Pb(col)).^2);
end


figure
plot(xCenters, xVariog, '-or')
hold on
plot(yCenters, yVariog, '-ob')
title('Directional variogram of P')
legend('E-W', 'N-S')



