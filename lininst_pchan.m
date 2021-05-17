
 % for jjj=[1 12 14 15 16 2 13 9], [~,casename] = system(['echo -n `basename ~/fms-backup/b',num2str(jjj),'[a-z]*`']); lininst_pchan; end
 % ${MiMAROOT}/src/shared/constants/constants.f90

%casename='b1ctrl';
if (exist('casename','var')==0), casename=getenv('CASENAME'); end
disp(casename);
bk = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'bk');
phalf = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'phalf');
bk_f = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'pfull') ./ phalf(end);
latf = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'lat');
latf = latf(end/2+1:end);
psavg = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'ps');
psavg = ( psavg(end/2+1:end)+psavg(end/2:-1:1) )/2;
tfull = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'temp');
ufull = ncread(['ensemble-wise/zyorg_',casename,'.nc'],'ucomp');
tfull = ( tfull(end/2:-1:1,:)+tfull(end/2+1:end,:) )/2;
ufull = ( ufull(end/2:-1:1,:)+ufull(end/2+1:end,:) )/2;
nlatf = numel(latf);
nArr = 0:15;

prm.drag = 1./86400*max(0,bk_f-0.7)/0.3;
prm.gravity = 9.80;
p_yp = psavg.*bk_f(:)';
theta = tfull.*(1e5./p_yp).^(2./7.);
strat = -287.04./(psavg.*bk(2:end-1)').*(psavg.*bk(2:end-1)'/1e5).^(2./7.).*diff(theta,1,2)./diff(p_yp,1,2);

kLininst = zeros(size(latf));
GrowthRate = zeros(size(latf));
evec = zeros([numel(bk_f) nlatf]);
Growth_yn = zeros([numel(nArr) nlatf]);
LengthArr = [0.1:0.1:4.5]*1e6;
for jj=1:numel(latf)
  beta = 2*7.292e-5*cos(latf(jj)*pi/180)/6371.0e3;
  f = 2*7.292e-5*sin(latf(jj)*pi/180);
%  [kLininst(jj),GrowthRate(jj)] = fminbnd(@(k) -real(lininst(beta,f,ufull(jj,:),psavg(jj)*bk_f,287.04*tfull(jj,end)/psavg(jj)/bk_f(end),strat(jj,:),psavg(jj)*bk,k/1e6,prm)) ,0.1/6371.0e3*1e6,100/6371.0e3*1e6);
  wrk(:,jj) = real(lininst(beta,f,ufull(jj,:),psavg(jj)*bk_f,287.04*tfull(jj,end)/psavg(jj)/bk_f(end),strat(jj,:),psavg(jj)*bk,2*pi./LengthArr,prm));  %alternative??
  [~,ind] = max(wrk(:,jj));
  [kLininst(jj),GrowthRate(jj)] = fminbnd(@(k) -real(lininst(beta,f,ufull(jj,:),psavg(jj)*bk_f,287.04*tfull(jj,end)/psavg(jj)/bk_f(end),strat(jj,:),psavg(jj)*bk,k/1e6,prm)) ,2*pi./LengthArr(ind)*0.9*1e6,2*pi./LengthArr(ind)*1.1*1e6);
  [~,evec(:,jj)] = lininst(beta,f,ufull(jj,:),psavg(jj)*bk_f,287.04*tfull(jj,end)/psavg(jj)/bk_f(end),strat(jj,:),psavg(jj)*bk,kLininst(jj)/1e6,prm);
  [Growth_yn(:,jj),~] = lininst(beta,f,ufull(jj,:),psavg(jj)*bk_f,287.04*tfull(jj,end)/psavg(jj)/bk_f(end),strat(jj,:),psavg(jj)*bk,nArr/6371.0e3/cos(latf(jj)*pi/180),prm);  %integer n
end
kLininst = kLininst/1e6;
GrowthRate = -GrowthRate;

%x{
addpath('/n/home05/pchan/bin');
figure;
set(gcf,'units','inches','position',[1 1 11 8.5], 'paperUnits','inches','papersize',[11 8.5],'paperposition',[0 0 11 8.5]);
%pcolorPH([0:0.05:1]*4,bk_f*1000,real(bb.*exp(1i*[0:0.05:1]*2*pi)));colorbar;grid on;axis ij;%caxis([1 1.001]);
subplot(2,1,1);
pcolorPH(-latf,bk_f*1000,real(evec));colorbar;grid on;axis ij;
colormap(gca,b2rCD(10)); caxis(0.4*[-1 1]);% caxis(max(max(abs(real(evec))))*[-1 1]);
subplot(2,1,2);
pcolorPH(-latf,bk_f*1000,imag(evec));colorbar;grid on;axis ij;
colormap(gca,b2rCD(10)); caxis(0.4*[-1 1]);% caxis(max(max(abs(real(evec))))*[-1 1]);

figure;
set(gcf,'units','inches','position',[1 1 6 4.5], 'paperUnits','inches','papersize',[6 4.5],'paperposition',[0 0 6 4.5]);
%pcolorPH(-latf,LengthArr,wrk./GrowthRate(:)');colorbar;grid on;caxis([1 1.001]);
wrkplot=(wrk./GrowthRate(:)'-0.9)*2000; wrkplot(wrkplot>200)=255;
image(-latf,LengthArr,wrkplot);axis xy;%colorbar;grid on;
hold on; plot(-latf,2*pi./kLininst,'m-o')
print(gcf, '-dpdf',['ensemble-wise/lininst_',casename,'.pdf']);
%x}
maxwrk = max(wrk)';
disp(latf( maxwrk>GrowthRate)')
disp(latf( maxwrk>GrowthRate*1.0001)')
%maxwrk( maxwrk>GrowthRate)

%disp('return'); return;
fn_savenc = ['ensemble-wise/lininst_',casename,'.nc'];
 system(['rm ',fn_savenc]);
nccreate(fn_savenc,'rGrowth_yn','Dimensions',{'kx',numel(nArr),'lat',nlatf},'DataType','single','Format','netcdf4')
nccreate(fn_savenc,'iGrowth_yn','Dimensions',{'kx',numel(nArr),'lat',nlatf},'DataType','single')
nccreate(fn_savenc,'kLininst','Dimensions',{'lat',nlatf},'DataType','single')
nccreate(fn_savenc,'lLininst','Dimensions',{'lat',nlatf},'DataType','single')
nccreate(fn_savenc,'GrowthRate','Dimensions',{'lat',nlatf},'DataType','single')
nccreate(fn_savenc,'lat','Dimensions',{'lat',nlatf},'DataType','single')
nccreate(fn_savenc,'kx','Dimensions',{'kx',numel(nArr)},'DataType','single')
%nccreate(fn_savenc,'p','Dimensions',{'p',npf},'DataType','single')
ncwriteatt(fn_savenc,'lat','units','degrees_N')
ncwriteatt(fn_savenc,'lat','long_name','latitude')
%ncwriteatt(fn_savenc,'p','units','hPa')
%ncwriteatt(fn_savenc,'p','long_name','pressure')
ncwrite(fn_savenc,'rGrowth_yn',real(Growth_yn(:,end:-1:1)))
ncwrite(fn_savenc,'iGrowth_yn',imag(Growth_yn(:,end:-1:1)))  % -ve for eastward
ncwrite(fn_savenc,'kLininst',kLininst(end:-1:1))
ncwrite(fn_savenc,'lLininst',2*pi./kLininst(end:-1:1))
ncwrite(fn_savenc,'GrowthRate',GrowthRate(end:-1:1))
ncwrite(fn_savenc,'lat',-latf(end:-1:1))
ncwrite(fn_savenc,'kx',nArr)
%ncwrite(fn_savenc,'p',presf)

% vim: set fdm=marker foldmarker=%{,%}:

