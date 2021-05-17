 % generate gridfrac.nc
 % export infile=d100/member/nc.2cospectrum_001.nc

 % n=0 excluded. Assumed n>0, then c=w/n monotonic with w. Integration of w is thus ok.
 % n<ds(3)/2+0.5 is upper limit.
 % w=0 retained. Ok for zonally sym model?
 % w>-ds(1)/2-0.5 is lower limit.
 % w<+ds(1)/2-0.5 is upper limit.

infile = getenv('infile');
ds = double(ncreadatt(infile,'/','ds')); ds(1)=[];
dt = double(ncreadatt(infile,'/','dt'));
Lx = double(ncreadatt(infile,'/','Lx'));
%ds(1)=10; ds(3)=16;
%ds(1)=100; ds(3)=192;
%ds(1)=400; ds(3)=128;
%dt=86400; Lx=6.371e6*2*pi;
%dt=86400/4; Lx=24e6;

%w_wn = repmat([0:ds(1)-1]',[1 ds(3)/2+1]);
w_wn = repmat(-ds(1)/2.+mod(ds(1)/2.+[0:ds(1)-1]',ds(1)),[1 ds(3)/2+1]);
n_wn = repmat([0:ds(3)/2],[ds(1) 1]);

c_bnd = [-20:55]';
phasespeed = (c_bnd(1:end-1)+c_bnd(2:end))/2;
%c_cwn = repmat(c_bnd(:),[1 ds(1) ds(3)/2+1]);
%w_cwn = repmat(-ds(1)/2.+mod(ds(1)/2.+[0:ds(1)-1],ds(1)),[numel(c_bnd) 1 ds(3)/2+1]);
%n_cwn = repmat(reshape(0:ds(3)/2,[1 1 ds(3)/2+1]),[numel(c_bnd) ds(1) 1]);

%frac_cwn = integral(@(nOffset) max(0,min(1, ...
%    c_cwn.*(n_cwn+nOffset)/Lx*dt*ds(1) -w_cwn+0.5)),-0.5,0.5,'ArrayValued',true);
cumsum_cwn = nan([numel(c_bnd) ds(1) ds(3)/2+1]);
for mc = 1:numel(c_bnd)
    cumsum_cwn(mc,:,:) = integral(@(nOffset) max(0,min(1, ...
        c_bnd(mc)*(n_wn+nOffset)/Lx*dt*ds(1) -w_wn+0.5)),-0.5,0.5,'ArrayValued',true);
    disp(mc);
end
cumsum_cwn(cumsum_cwn>1)=1;
%cumsum_cwn(cumsum_cwn<1e-7)=0;
cumsum_cwn(:,:,1)=0;  % n=0 excluded

frac_cwn = diff(cumsum_cwn,1,1);
frac_cwn(frac_cwn<0) = 0;

fn_savenc = ['gridfrac.nc'];
 system(['rm ',fn_savenc]);
nccreate(fn_savenc,'frac_cwn','Dimensions',{'kx',ds(3)/2+1,'freq',ds(1),'phasespeed',numel(phasespeed)},'DataType','single','Format','netcdf4')
nccreate(fn_savenc,'phasespeed','Dimensions',{'phasespeed',numel(phasespeed)},'DataType','single')
%ncwriteatt(fn_savenc,'lat','units','degrees_N')
%ncwriteatt(fn_savenc,'lat','long_name','latitude')
ncwrite(fn_savenc,'frac_cwn',permute(frac_cwn,3:-1:1))
ncwrite(fn_savenc,'phasespeed',phasespeed)


%{
%f_c_wn = @(jw,in)( -ds(1)/2.+mod(ds(1)/2.+jw,ds(1)) )*2*pi/ds(1)/dt ./in*6.371e6;
%f_c_wn = @(jw,in) jw*2*pi/ds(1)/dt ./in*6.371e6;
%f_c_wn = @(wOffset,nOffset) (w_wn+wOffset)*2*pi/ds(1)/dt ./(n_wn+nOffset)*6.371e6;
f_c_wn = @(wOffset,nOffset) repmat( reshape((w_wn+wOffset)*2*pi/ds(1)/dt ./(n_wn+nOffset)*6.371e6 ...
    ,[1 ds(1) ds(3)/2+1]),[numel(c_bnd) 1 1])<c_cwn;
w1_cwn = repmat(c_bnd,[1 ds(1) ds(3)/2+1]).*repmat(n_wn-0.5,;

code_cwn = 20*f_c_wn(-0.5,-0.5)+10*f_c_wn(0.5,-0.5)+2*f_c_wn(-0.5,0.5)+1*f_c_wn(0.5,0.5);
code_cwn(:,:,1)=0;

%unique(code_cwn), unique(ismember(code_cwn,[0 20 2 22 3 30 23 32 33])),
%}

% vim: set fdm=marker foldmarker=%{,%}:

