clear
close all

%%
d = importdata('obs_u.mat');
G = importdata('obs_G.mat');
L = importdata('laplacian.mat');


d_ = [d;zeros(320,1)];

n = 50;
alpha_vec = exp(linspace(-5,2,n));
misfit_vec = zeros(1,n);
%roughness_vec = zeros(1,n);
for ii = 1:n 
    G_ = [G;alpha_vec(ii)*L];
    m = (G_'*G_)\(G_'*d_);
    
    v = d - G*m;
    misfit_vec(ii) = v'*v;
    
    %p = L*m;
    %roughness_vec(ii) = sum(p)/640;
end

figure;
plot(log(alpha_vec),misfit_vec,'k.-');

% figure;
% plot(roughness_vec,misfit_vec,'k.-');

%m = fminunc(@(x) (d-G*x)'*(d-G*x) + 100*(L*x)'*(L*x), 1*randn(320,1));

%save('ruptures','-v7.3','M');
%save('data','-v7.3','M');

%%
G_ = [G;6*L];
m = (G_'*G_)\(G_'*d_);
figure;
imagesc(reshape(m,8,40));
colorbar;
colormap('jet');
%caxis([-10,20]);
wysiwyg(gcf,20,3);

G_ = [G;0.3*L];
m = (G_'*G_)\(G_'*d_);
figure;
imagesc(reshape(m,8,40));
colorbar;
colormap('jet');
%caxis([-10,20]);
wysiwyg(gcf,20,3);

G_ = [G;1e-8*L];
m = (G_'*G_)\(G_'*d_);
figure;
imagesc(reshape(m,8,40));
colorbar;
colormap('jet');
%caxis([-10,20]);
wysiwyg(gcf,20,3);

% more noise in data AND possible noise in slip distribution!