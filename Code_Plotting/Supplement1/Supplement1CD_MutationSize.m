MutationPDF = [3.1, 1.2, 3.08];
X = normrnd(MutationPDF(1), MutationPDF(2), 1, 100000);
dE = log10(exp(1))*(exp(X) - MutationPDF(3));
a = sum(10.^(-dE)>10);
b = sum(10.^(-dE)>1);

figure
histogram(-dE)
x = linspace(-4, 6, 10001);
y = pdf('normal', x, MutationPDF(1), MutationPDF(2));
xprime = log10(exp(1))*(exp(x)-MutationPDF(3));
yprime = y./(log10(exp(1))*exp(x));

figure
t = tiledlayout(1,2,'padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 10, 5];
nexttile;
plot(-xprime, yprime, 'color','k')
xlim([-3, 2])
xlabel('Change in affinity (-log_{10}K_d)', 'fontsize', 8)
ylabel('Probability Density', 'fontsize', 8)
grid on
xticks([-8:2])
nexttile;
cdf = cumtrapz(xprime, yprime);
plot(-xprime, 1-cdf, 'color', 'k')
xlim([-2, 2])
xlabel('Change in affinity (-log_{10}K_d)', 'fontsize', 8)
ylabel('Cumulative Probability', 'fontsize', 8)
grid on
f = gcf;
exportgraphics(f,fullfile('..','figures','Supplement1c_MutationScheme.pdf'),'ContentType','vector',...
            'BackgroundColor','none');
%%
MutationPDF = [3.1, 1.2, 3.08];
rhos = [0.4:0.02:0.8, 0.905:0.005:1];
double_ben = zeros(size(rhos));
N = 10^7;
ben1 = zeros(size(rhos));
ben2 = zeros(size(rhos));

m = 5;
for i=1:m
for rho = rhos
    sigma = [1,rho;rho,1];
    X = MutationPDF(1) + mvnrnd([0,0],MutationPDF(2)^2*sigma, N);
    dE = log10(exp(1))*(exp(X) - MutationPDF(3));
    double_ben(rho==rhos)=double_ben(rho==rhos)+sum(dE(:,1)<0 & dE(:,2)<0);
    ben1(rho==rhos)=ben1(rho==rhos)+sum(dE(:,1)<0);
    ben2(rho==rhos)=ben2(rho==rhos)+sum(dE(:,2)<0);
end
end
double_ben_frac = double_ben/(m*N);
ben1frac = ben1/(m*N);
ben2farc = ben2/(m*N);
%%
figure
t = tiledlayout(1,1,'padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3,3,5.5,5];

nexttile;
plot(rhos,double_ben_frac./ben1frac, 'color', 'k');
xticks(0.4:0.2:1)
xlabel('\rho', 'fontsize', 8);
ylabel({'Fraction of beneficial mutations that', 'are beneficial for binding both strains'}, 'fontsize', 8)
grid on
f = gcf;
savefig(f, fullfile('..','figures','Supplement1d_double_beneficial.fig'))
exportgraphics(f,fullfile('..','figures','Supplement1d_double_beneficial.pdf'),'ContentType','vector',...
            'BackgroundColor','none');