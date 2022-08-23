dose3clones=[0.300946413, 1.833449542, 12.70520533, 937.1275442, 1000, 1000, 0.565667551, 13.28638369, 15.51126244, 19.97886876, 28.20299352, 34.7263686, 46.86562478, 79.7205409, 120.2399052, 717.2925921, 866.0502971, 1000, 1000, 1000, 1000, 1000, 1000, 1000];
dose3egc=[0.926815928,522.5195881	7.443800546	300.5836239	270.9899404	44.66065782	36.37397655	931.3176355];
dose3sing=[20.27258171, 64.07420332, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 2.293360808, 121.3005791];

mo5clones=[26.42266637, 28.3591188, 49.80031108, 104.730714, 105.3683844, 124.224531, 159.3412508, 666.6414773, 681.2685339, 1000, 1000, 1000, 1000, 1000, 1000, 1000,522.5195881];
mo5egc= [896.6985098,948.5173907, 694.8,931.3176355];
mo5sing=[1	1000	1000	1000	126.1521802	1000	493.9062623	1000	1000	8.482500466	418.7675448	1000	1000	63.49200699	1000	2.65763454	108.2521313	1000	1000	1000	1000	1000	1000	1000	1000	1000];

mo13egc=[1000 1000 918.8582911 920.379299 1000];
mo13rando=[1000	1000	1000	1000	1000	1000	1000	1000	1000	1000	160.2916901	1000	1000	1000	1000	1000	645.9524383	1000	1000	392.1645977];


[f,x]=ecdf(-log([mo5egc, mo13egc])/log(10)); [f2,x2]=ecdf(-log(mo5clones)/log(10));
[f3,x3]=ecdf(-log(mo5sing)/log(10)); [f4,x4]=ecdf(-log(dose3egc)/log(10));
[f5,x5]=ecdf(-log(dose3clones)/log(10)); [f6,x6]=ecdf(-log(dose3sing)/log(10));
[f7,x7]=ecdf(-log(mo13egc)/log(10)); [f8,x8]=ecdf(-log(mo13rando)/log(10));



 xlow=.07;  ylow=.17;  xwidth=.35;  yheight=.60;
 figure('Position', [10 10 1000 400]*.70) 
subplot('Position',[xlow ylow xwidth yheight])

% figure
h=plot(x,f,x2,f2,x3,f3,'Linewidth',2);
xlabel('-log_1_0 (IC_5_0 Omicron)','Fontweight','Normal')
ylabel('Cumulative Fraction','Fontweight','Normal')
 h(3).Color=sscanf('95b5db' ,'%2x%2x%2x',[1 3])/255;
      h(2).Color=[0, 0.4470, 0.7410];
      h(1).Color=[0.8500, 0.3250, 0.0980];
line([mean(-log([mo5egc, mo13egc])/log(10)), mean(-log([mo5egc, mo13egc])/log(10))],[0,1],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',1.5)
line([mean(-log(mo5clones)/log(10)), mean(-log(mo5clones)/log(10))],[0,1],'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',1.5)
line([mean(-log(mo5sing)/log(10)), mean(-log(mo5sing)/log(10))],[0,1],'Color',sscanf('95b5db' ,'%2x%2x%2x',[1 3])/255,'LineStyle','--','LineWidth',1.5)
legend('EGC Clones','Other 5 month Clones','5 month singlets', 'location', 'east')
xlim([-3.1,1])

 t=title('A');  t.FontSize = 12;
 set(t,'position',get(t,'position')-[2.5 0 0])
 set(gca, 'FontName', 'Arial')
 subplot('Position',[.5+xlow ylow xwidth yheight])

 
 
 
% figure
h=plot(x4,f4,x5,f5,x6,f6,'Linewidth',2);
xlabel('-log_1_0 (IC_5_0 Omicron)','Fontweight','Normal')
ylabel('Cumulative Fraction','Fontweight','Normal')
       h(3).Color=sscanf('95b5db' ,'%2x%2x%2x',[1 3])/255;
      h(2).Color=[0, 0.4470, 0.7410];
      h(1).Color=[0.8500, 0.3250, 0.0980];
line([mean(-log(dose3egc)/log(10)), mean(-log(dose3egc)/log(10))],[0,1],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',1.5)
line([mean(-log(dose3clones)/log(10)), mean(-log(dose3clones)/log(10))],[0,1],'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',1.5)
line([mean(-log(dose3sing)/log(10)), mean(-log(dose3sing)/log(10))],[0,1],'Color',sscanf('95b5db' ,'%2x%2x%2x',[1 3])/255,'LineStyle','--','LineWidth',1.5)
legend('EGC Clones','Other Vax3 Clones' ,'Vax3 singlets', 'location', 'east')
xlim([-3.1,1])
 t=title('B');  t.FontSize = 12;
 set(t,'position',get(t,'position')-[2.5 0 0])
 set(gca, 'FontName', 'Arial')
 
exportgraphics(gcf, fullfile('..','figures','Figure4.pdf'), 'Resolution', 600)
 
 
 %1.3 month figure
 
 
figure
h=plot(x7,f7,x8,f8,'Linewidth',2);
xlabel('-log_1_0 (IC_5_0 Omicron)','Fontweight','Normal')
ylabel('Cumulative Fraction','Fontweight','Normal')
      h(2).Color=[0, 0.4470, 0.7410];
      h(1).Color=[0.8500, 0.3250, 0.0980];
line([mean(-log(mo13egc)/log(10)), mean(-log(mo13egc)/log(10))],[0,1],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',1.5)
line([mean(-log(mo13rando)/log(10)), mean(-log(mo13rando)/log(10))],[0,1],'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',1.5)
legend('EGC Clones','Other 1.3 mo B cells', 'location', 'east')
xlim([-3.1,1])











%%Leerang can ignore this, it's just stuff I had for stats.
x1=mean(log(dose3egc));
s1=std(log(dose3egc));
x2=mean(log([mo5egc, mo13egc]));
s2=std(log([mo5egc, mo13egc]));
x3=mean(log([mo5sing, mo5clones]));
s3=std(log([mo5sing, mo5clones]));
% x2=mean(log(dose3clones));
% s2=std(log(dose3clones));
t=(x1-x2)/sqrt(s1^2 /length(dose3egc) + s2^2 /length([mo5egc, mo13egc]));
t2=(x1-x3)/sqrt(s1^2 /length(dose3egc) + s3^2 /length([mo5sing, mo5clones]));
t3=(x2-x3)/sqrt(s2^2 /length([mo5egc, mo13egc]) + s3^2 /length([mo5sing, mo5clones]));
% t=(x1-x2)/sqrt(s1^2 /length(dose3egc) + s2^2 /length(dose3clones));

exp(mean(log([mo5egc, mo13egc])))
exp(mean(log(dose3egc)))


exp(mean(log([mo5clones,mo5sing])))
exp(mean(log([dose3clones,dose3sing])))
