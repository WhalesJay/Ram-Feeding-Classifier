function [gapearea,ci_gapes,width,ci_width,lnth] = getgape(age,lnth) 
% get estimated gape area for NARW of a given age
% inputs: 
%      age: age of animal (years)
%      lnth: length (cm) if known -- default is age-based from Moore et al. 2005
% outputs:
%      gape: area in m^2
%      ci_gape: 95% confidence interval around gape estimate, propagated
%      width: width of snout based on Miller et al
%      ci_width: 95% confidence interval around width estimate
%      lnth: length (either estimated or output from Moore)

if nargin <2
% estimate length of whale from age (Moore et al. 2004)
[lnth,ci_length] = FortuneAgeLength(age);
else
ci_length = [lnth lnth]; 
end

% load length/witdh data from Carolyn Miller: Eubalaena glacialis and
% Eubalaena australis calves and moms
load('EgEaWidthsCAMthesis')
% fit 10% width to length 
wd10 = fit(data(:,1),data(:,2),'poly1'); 
% estimate witdh from length
width = feval(wd10,lnth); 
ci_width = predint(wd10,lnth,0.95,'functional','off'); % confidence interval 

%% add  baleen length
%Best and Schell 1996, table 1
bodylength = [13.1 10.06 14.3 9.25 4.6 13.76 14.8 10.31 11.76 11.23 10.66];
platelength = [1.75 1.025 1.90 1.05 .12 2.375 2.025 1.125 .98 .82 1.36]; % (m)

% Omura (1958), Klumov (1962), Matthews (1938) as reported in Omura 1969
Antarctic = [6.4664073	0.1593402; 13.501248	1.8153068; 14.4106	1.9719481]; 
NPRW = [10.777474	0.39463958
11.657088	0.89067054
12.418416	1.9719481
12.612979	0.98490834
14.112836	2.0025098
14.711208	2.0512502
15.111184	1.8767136
15.203412	2.0911121
15.426739	2.1613517
16.039192	2.3819911
16.154743	2.3812551
16.309784	2.0458448
16.621752	2.1919591
16.622633	2.5980382
16.416725	2.8286684
17.025084	2.519035
17.084837	2.595094
17.115084	2.6474535
17.441118	2.5928245
17.820353	2.3610904
18.331057	1.9995272]; 

% figure(21), clf, hold on 
% plot(bodylength,platelength,'v')
% plot(NPRW(:,1),NPRW(:,2),'^')
% plot(Antarctic(:,1),Antarctic(:,2),'o')
% 
% xlabel('Body Length (m)')
% ylabel('Longest Baleen Plate Length (m)')
% adjustfigurefont 

allBlength = vertcat(bodylength',NPRW(:,1),Antarctic(:,1));
allPlength = vertcat(platelength',NPRW(:,2),Antarctic(:,2));

% fit the model 
bln = fit(allBlength,allPlength,'poly1'); 
% estimate baleen length from body length
baleen = feval(bln,lnth/100);
ci_baleen = predint(bln,lnth/100,0.95,'functional','off');

errorProp % propagate errors into gape 

gapearea = (width*baleen)/2; % m^2

ci_gapes = del_all*gapearea; 
