function [length_out,ci_length] = FortuneAgeLength(age_in)

% obtain fit and prediction interval for published relationship in Fortune et
% al. 2005 Morphometry and gross morphology right whales... 

% data from Table 4 Moore et al 2005 Morphometry
lnth = [996 1048 1093 1132 1167 1194 1218 1239 1256 1272 1285 1296 1306 1315 1322 1328 1333 1338 1342 1345 1348 1351 1353 1355 1357 1358 1359 1360 1361 1362 1030];
age = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];

% find Nans in age and remove 
ii = isnan(age); 
lnth(ii) = []; 
age(ii) = []; 
% remove age zero
ii = find(age ~= 0);


% plot
% figure(1), clf, hold on 
% plot(age,lnth,'o')
% plot Michael's Fit
% whaleAge = 1:30;
% lnthFit(whaleAge) = 1011.033+320.501*log10(whaleAge); % MOORE ET AL 2004
% plot(whaleAge,lnthFit)

% fit
ft = fit(log10(age(ii)'),lnth(ii)','poly1');
% good, same fit as Michael's 

% prediction interval 
length_out = feval(ft,log10(age_in)); 
ci_length = predint(ft,log10(age_in),0.68,'functional','off');

return 
% for paper 
paper_lengths_out = feval(ft,log10([4 3 18 4 11 2 11 3 3 4]));
paper_ci_length = predint(ft,log10([4 3 18 4 11 2 11 3 3 4]),0.68,'functional','off');
sds = (paper_ci_length(:,1)-paper_lengths_out); % standard deviations 
figure(99), clf, hold on
errorbar([4 3 18 4 11 2 11 3 3 4],paper_lengths_out/100,sds/100,'o') % plot in m
plot(age,lnth/100,'o')
xlabel('Age'), ylabel('Body Length (m)') 

% plot(whaleAge,ci_length(:,1))
% plot(whaleAge,ci_length(:,2))
% xlabel('Age'), ylabel('Length (cm)'), adjustfigurefont