clear all
load('DataHW2.mat')


%% Problem 1
% (3 points)

%%
% a) (0.5 point)
% We can see the famous football player Lionel Messi.

figure(1)
imagesc(M)
colormap('gray')
colorbar()
title('Lionel Messi')

%%
% b) (1.5 point)
% We choose a moving window of size 30x30. We use the function nlfilter 
% which pads the data points beyond the edges with 0s. We can see that this
% solution creates some artifacts at the border of the image and is not
% ideal. However we have enough data here so we have good hope that it
% won't affect too much the analysis. 
% 

figure(2)
M_mean = nlfilter(M, [30 30], @(x) mean(x(:)));
imagesc(M_mean)
colormap('gray')
colorbar()
title('Moving-window mean')

figure(3)
M_std = nlfilter(M, [30 30], @(x) std(x(:)));
imagesc(M_std)
colormap('gray')
colorbar()
title('Moving-window std')

%%
% c) (1 point)
% According to the correlation coefficient there is a significant positive
% correlation between the mean and the standard deviation. 
% Overall we can conclude that there is a proportional effect. 

figure(4)
scatter(M_mean(:), M_std(:), .01, 'b.')
xlabel('mean')
ylabel('std')
hold on

p = polyfit(M_mean(:), M_std(:),1);
x = -6:250;
y = polyval(p, x);
plot(x, y, 'r', 'LineWidth', 2)
hold off

r = corr(log(M_mean(:)+1), log(M_std(:)+1))

% This can be also seen from the scatter plot, even though there is a very 
% important noise.

figure(5)
[M_mean_sort, ind] = sort(M_mean(:));
M_std_sort = M_std(:);
M_std_sort = M_std_sort(ind);
M_std_sort = reshape(M_std_sort, [], 10);
M_mean_sort = reshape(M_mean_sort, [], 10);
M_mean_bins = mean(M_mean_sort);
boxplot(M_std_sort, 'labels', M_mean_bins);

% The best way to see it however is by using a log transformation

figure(6)
scatter(log(M_mean(:)), log(M_std(:)), .01, 'b.')
xlabel('mean')
ylabel('std')
hold on

p = polyfit(log(M_mean(:)+.001), log(M_std(:)+.001),1);
x = -6:6;
y = polyval(p, x);
plot(x, y, 'r', 'LineWidth', 2)
hold off

%% Problem 2
% (5 points)
clear all
load('Data_HW1.mat')


%%
% a) (1 point)

figure(6)
Rmed = median(R(:));
Rind = R > Rmed;
imagesc(Rind)                               
cmap = jet(2);                       
colormap(cmap)
hold on
L = line(ones(2),ones(2), 'LineWidth',2);     
set(L,{'color'},mat2cell(cmap,ones(1,2),3));  
legend('<= median(R)','> median(R)')  
hold off

figure(7)
Tmed = median(T(:));
Tind = T > Tmed;
imagesc(Tind)                            
cmap = jet(2);                       
colormap(cmap)
hold on
L = line(ones(2),ones(2), 'LineWidth',2);     
set(L,{'color'},mat2cell(cmap,ones(1,2),3));  
legend('<= median(T)','> median(T)')  
hold off

%%
% b) 

figure(8)
[C,h] = contour(T);
clabel(C, h)
title('Contour map of T')

%% 
% c)

figure(9)

h_scatter = @(n) scatter(reshape(R(1:end-n, 1:end-n), 1, []), ...
                         reshape(R(n+1:end, n+1:end), 1, []));

for i=1:6
    subplot(2,3,i); h_scatter(i); title(sprintf('(%d,%d)-scatterplot',[i,i])); 
    ylim([0 70]); xlabel('R(t)') ; ylabel('R(t+h)')
end

%%
% d)
n = size(R, 1)-1;
m = size(R, 2)-1;

%%% correlogram %%%
h_correlogram = @(i, j) corr(reshape(R(1:end-i, 1:end-j), [], 1),...
                             reshape(R(i+1:end, j+1:end), [], 1));
correlogram = zeros([n, m]);          
for i=1:n
    for j=1:m
        correlogram(i, j) = h_correlogram(i, j);    
    end
end

figure(10)
imagesc(correlogram)
colormap('cool')
colorbar()
xlabel('x')
ylabel('y')
title('Correlogram')

%%% covariogram %%%
h_covariogram = @(i, j) cov(reshape(R(1:end-i, 1:end-j), [], 1),...
                            reshape(R(i+1:end, j+1:end), [], 1));
covariogram = zeros([n, m]);          
for i=1:n
    for j=1:m
        h_cov = h_covariogram(i, j);
        covariogram(i, j) = h_cov(1, 2);
    end
end

figure(11)
imagesc(covariogram)
colormap('cool')
colorbar()
xlabel('x')
ylabel('y')
title('Covariogram')

%%

%%% variogram %%%

variogram = zeros([n, m]);          
for i=1:n
    for j=1:m
        diff = reshape(R(1:end-i, 1:end-j), [], 1) ...
               - reshape(R(i+1:end, j+1:end), [], 1);
        variog = .5*mean(diff.^2);
        variogram(i, j) = variog;
    end
end

figure(12)
imagesc(variogram)
colormap('cool')
colorbar()
xlabel('x')
ylabel('y')
title('Variogram')


%% 
% e)
% We observe that in the (h,-h) direction, the data set remains much more
% corelated than in the (h,h) direction. This is coherent with what 
% indicates the correlogram. Plotting the R dataset, we can clearly
% distinguish structures in the (h,-h) direction.

figure(13)

h_scatter = @(n) scatter(reshape(R(1:end-n, n+1:end), 1, []), ...
                         reshape(R(n+1:end, 1:end-n), 1, []));

for i=1:6
    subplot(2,3,i); h_scatter(i); title(sprintf('(%d,-%d)-scatterplot', i,i)); 
    ylim([0 70]); xlabel('R(t)') ; ylabel('R(t+h)')
end

figure(14)
imagesc(R)
title('R values')


%% 
% d)
% Case (h, h)
% sill: 98
% range : 30
% nugget : 1

N = 50;
variogram = zeros([N, 1]);          
for i=1:N
    diff = reshape(R(1:end-i, 1:end-i), [], 1) ...
           - reshape(R(i+1:end, i+1:end), [], 1);
    variog = .5*mean(diff.^2);
    variogram(i) = variog;
end

figure(14)
plot(variogram)

%%
% Case (h, -h)
% sill : ?
% range : ?
% nugget : 0

N = 60;
variogram = zeros([N, 1]);          
for i=1:N
    diff = reshape(R(1:end-i, 1+i:end), [], 1) ...
           - reshape(R(i+1:end, 1:end-i), [], 1);
    variog = .5*mean(diff.^2);
    variogram(i) = variog;
end

figure(15)
plot(variogram)

%% Problem 3
% (2 points)

clear all
load('Data_HW1.mat')

%%
% (a)
figure(16)

h_scatter = @(n) scatter(reshape(R(1:end-n, 1:end-n), 1, []), ...
                         reshape(T(n+1:end, n+1:end), 1, []));
                     

subplot(221); h_scatter(0); title('cross (0,0)-scatterplot'); 
xlabel('R(t)') ; ylabel('T(t+h)')
subplot(222); h_scatter(1); title('cross (1,1)-scatterplot'); 
xlabel('R(t)') ; ylabel('T(t+h)')
subplot(223); h_scatter(2); title('cross (2,2)-scatterplot'); 
xlabel('R(t)') ; ylabel('T(t+h)')
subplot(224); h_scatter(3); title('cross (3,3)-scatterplot'); 
xlabel('R(t)') ; ylabel('T(t+h)')



h_corr = @(n) corr(reshape(R(1:end-n, 1:end-n),[],1), ...
                         reshape(T(n+1:end, n+1:end),[],1));
                     
corr00 = h_corr(0)
corr11 = h_corr(1)
corr22 = h_corr(2)
corr33 = h_corr(3)

%%
% (b)

m = size(R, 2)-1;

%%%%%%%%%%%%% along (h,h) axis %%%%%%%%%%%%%%

%%% correlogram %%%
h_correlogram = @(h) corr(reshape(T(1:end-h, 1:end-h), [], 1),...
                          reshape(R(h+1:end, h+1:end), [], 1));
correlogram = zeros([m,1]);          
for h=1:m
    correlogram(h) = h_correlogram(h);    
end

figure(17)
plot(correlogram)
xlabel('h')
ylabel('correlation')
title('Correlogram')

%%% variogram %%%
variogram = zeros([m,1]);          
for h=1:m
    diff1 = reshape(R(1:end-h, 1:end-h), [], 1) ...
          - reshape(R(h+1:end, h+1:end), [], 1);
    diff2 = reshape(T(1:end-h, 1:end-h), [], 1) ...
          - reshape(T(h+1:end, h+1:end), [], 1);
    variog = .5*mean(diff1 .* diff2);
    variogram(h) = variog;   
end

figure(18)
plot(variogram)
xlabel('h')
ylabel('variance')
title('Variogram')


%%
%%%%%%%%%%%%% along (h, -h) axis %%%%%%%%%%%

%%% correlogram %%%
h_correlogram = @(h) corr(reshape(T(1:end-h, h+1:end), [], 1),...
                          reshape(R(h+1:end, 1:end-h), [], 1));
correlogram = zeros([m,1]);          
for h=1:m
    correlogram(h) = h_correlogram(h);    
end

figure(19)
plot(correlogram)
xlabel('h')
ylabel('correlation')
title('Correlogram')

%%% variogram %%%
variogram = zeros([m,1]);          
for h=1:m
    diff1 = reshape(R(1:end-h, h+1:end), [], 1) ...
          - reshape(R(h+1:end, 1:end-h), [], 1);
    diff2 = reshape(T(1:end-h, h+1:end), [], 1) ...
          - reshape(T(h+1:end, 1:end-h), [], 1);
    variog = .5*mean(diff1 .* diff2);
    variogram(h) = variog;   
end

figure(20)
plot(variogram)
xlabel('h')
ylabel('variance')
title('Variogram')


