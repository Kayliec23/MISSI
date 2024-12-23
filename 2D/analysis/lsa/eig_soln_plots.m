addpath /Users/Kayliec23/Desktop/plottingKit/ 
addpath /Users/Kayliec23/Desktop/plottingKit/addaxis6
addpath /Users/Kayliec23/Desktop/analysis
addpath /Users/Kayliec23/Desktop/MISSI/analysis/DiffusiveInstability-main

marg_h = [];
marg_w = [];
gap = [0.05 0.1];
fontsize = 24;

N2 = 1.8e-4; %%% Buoyancy frequency
[gr,lamY,lamZ,Ri0,ll,mm] = analytical_soln(N2,0,1);

figure;
subtightplot(1,3,1,gap)
[MM,LL] = meshgrid(mm,ll);
pcolor(LL,MM,gr);
shading flat;
cb = colorbar;
cb.Location = 'northoutside';
cb.LineWidth = 2;
caxis([-1 1]*1e-4);
cmocean('curl');
xticks([-0.099 0 0.099]); xticklabels([-0.1 0 0.1])
yticks([-0.99 0 0.99]); yticklabels([-1 0 1])
set(gca,'box','on','linewidth',2);
set(gca,'fontsize',fontsize);
xlabel('$l$','interpreter','latex');
ylabel('$m$','interpreter','latex');
pbaspect([1 1 1])

%%
N2vals = [1.3:0.025:3]*1e-4;
maxgr = 0*N2vals;
Ri0vals = 0*N2vals;
lamYvals = 0*N2vals;
lamZvals = 0*N2vals;
for n=1:length(N2vals)
  [gr,lamY,lamZ,Ri0,ll,mm] = analytical_soln(N2vals(n),0,0);
  maxgr(n) = max(gr(:));
  if (maxgr(n)>0)
    [i,j] = find(gr==maxgr(n));
    lamYvals(n) = lamY(i);
    lamZvals(n) = lamZ(j);    
  end
  Ri0vals(n) = Ri0;
end

lamYvals(maxgr<=0) = NaN;
lamZvals(maxgr<=0) = NaN;

% this needs some manual adjustments... not sure how addaxis works
sb = subtightplot(1,3,2,gap); hold on;
plot(Ri0vals,lamYvals,'linewidth',2);
ax = sb;
ax.YColor = [0.47,0.67,0.19];

addaxis(Ri0vals,maxgr,[0 12]*1e-5,'linewidth',2)
addaxis(Ri0vals,lamZvals,[0 75],'linewidth',2)

addaxislabel(1,'\lambda_{horizontal}^{max} (m)')
addaxislabel(2,'\sigma^{max} (s^-^1)')
addaxislabel(3,'\lambda_{vertical}^{max} (m)')
xlabel('Ri')
set(gca,'box','on','linewidth',2);
set(gca,'fontsize',fontsize);
pbaspect([1 1 1])


%%
N2vals = [0.5:5e-4:4]*10e-3; %[1e-6:3e-7:1.2e-4,1.2e-4:3e-6:4e-3]; %[1.3:0.025:3]*1e-4;
svals = [0.05e-2 : 0.1e-2 : 5e-2];
[SS,NN2] = meshgrid(svals,N2vals);

for n = 1:length(N2vals)
    for k = 1:length(svals)
        [gr,lamY,lamZ,Ri0,ll,mm] = analytical_soln(N2vals(n),svals(k),0);
        maxgrvals(n,k) = max(gr(:));
    end
end

%figure;
subtightplot(1,3,3,gap)
pcolor(NN2,SS,maxgrvals);
shading interp;
cb = colorbar;
cb.Location = 'northoutside';
cb.LineWidth = 2;
xlabel('N^2 (s^-^2)')
ylabel('slope')
set(gca,'box','on','linewidth',2);
set(gca,'fontsize',fontsize);
pbaspect([1 1 1])
