addpath /Users/Kayliec23/Desktop/plottingKit/ 
addpath /Users/Kayliec23/Desktop/plottingKit/addaxis6
addpath /Users/Kayliec23/Desktop/analysis
addpath /Users/Kayliec23/Desktop/MISSI/analysis/DiffusiveInstability-main


fontsize = 24;

N2 = 2e-4; %%% Buoyancy frequency
[gr,lamY,lamZ,Ri0] = analytical_soln(N2,0);

figure;
subplot(3,1,1)
[MM,LL] = meshgrid(lamZ,lamY);
pcolor(LL,MM,gr);
shading flat;
cb = colorbar;
cb.Location = 'northoutside';
caxis([-5 5]*1e-5);
cmocean('curl');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'box','on','linewidth',2);
set(gca,'fontsize',fontsize);
xlabel('\lambda_{horizontal} (m)');
ylabel('\lambda_{vertical} (m)');


N2vals = [1.3:0.025:3]*1e-4;
maxgr = 0*N2vals;
Ri0vals = 0*N2vals;
lamYvals = 0*N2vals;
lamZvals = 0*N2vals;
for n=1:length(N2vals)
  [gr,lamY,lamZ,Ri0] = analytical_soln(N2vals(n),0);
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
sb = subplot(3,1,2); hold on;
plot(Ri0vals,lamYvals,'linewidth',2);
ax = sb;
ax.YColor = [0.47,0.67,0.19];

addaxis(Ri0vals,maxgr,[0 12]*1e-5,'linewidth',2)
addaxis(Ri0vals,lamZvals,[0 75],'linewidth',2)

addaxislabel(1,'\lambda_{horizontal}^{max} (m)')
addaxislabel(2,'\sigma^{max} (s^-^1)')
addaxislabel(3,'\lambda_{vertical}^{max} (m)')




%%
N2vals = [1e-6:3e-7:1.2e-4,1.2e-4:3e-5:4e-3]; %[1.3:0.025:3]*1e-4;
svals = [0.05e-2 : 0.1e-2 : 5e-2];
[SS,NN2] = meshgrid(svals,N2vals);

for n = 1:length(N2vals)
    for k = 1:length(svals)
        [gr,lamY,lamZ,Ri0] = analytical_soln(N2vals(n),svals(k));
        maxgrvals(n,k) = max(gr(:));
%         if (maxgr(n,k)>0)
%             [i,j] = find(gr==maxgr(n,k));
%             lamYvals(n,k) = lamY(i);
%             lamZvals(n,k) = lamZ(j);
%         else
%             lamYvals(n,k) = NaN;
%             lamZvals(n,k) = NaN;
%             SS(n,k) = NaN;
%             NN2(n,k) = NaN;
%         end
    end
end

%figure;
subplot(3,1,3)
pcolor(NN2,SS,maxgrvals);
shading interp;
colorbar;
xlabel('N^2 (s^-^2)')
ylabel('slope')



