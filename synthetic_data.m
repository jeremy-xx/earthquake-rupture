clear
close all

%% fault
l = 40;
w = 8;
s = 110;

n_h = 40;
n_v = 8;
n_a = n_h*n_v;

l_n = l/n_h;
w_n = w/(n_v-1);

x1 = -sind(s)*l/2;
y1 = -cosd(s)*l/2;
x2 = sind(s)*l/2;
y2 = cosd(s)*l/2;
z1 = 0;
z2 = -w;

x = linspace(x1,x2,n_h);
y = linspace(y1,y2,n_h);
z = linspace(z1,z2,n_v);

X = zeros(1,n_a);
Y = zeros(1,n_a);
Z = zeros(1,n_a);

faults_X = zeros(4,n_a);
faults_Y = zeros(4,n_a);
faults_Z = zeros(4,n_a);

for ii = 1:n_h
    for jj = 1:n_v
        kk = (ii-1)*n_v + jj;
        
        X(kk) = x(ii);
        Y(kk) = y(ii);
        Z(kk) = z(jj);
        
        faults_X(1,kk) = x(ii) - sind(s)*l_n/2;
        faults_X(2,kk) = x(ii) + sind(s)*l_n/2;
        faults_X(3,kk) = faults_X(2,kk);
        faults_X(4,kk) = faults_X(1,kk);
        
        faults_Y(1,kk) = y(ii) - cosd(s)*l_n/2;
        faults_Y(2,kk) = y(ii) + cosd(s)*l_n/2;
        faults_Y(3,kk) = faults_Y(2,kk);
        faults_Y(4,kk) = faults_Y(1,kk);
        
        faults_Z(1,kk) = z(jj);
        faults_Z(2,kk) = faults_Z(1,kk);
        faults_Z(3,kk) = z(jj) - w_n;
        faults_Z(4,kk) = faults_Z(3,kk);
    end
end

% fault ruptures
r = linspace(-l/2,l/2,n_h);
S = 10*gaussian(z',r,0,-12,12) + 11*gaussian(z',0.5*r,2,3,15);
S = S(:)';

S = S + randn(1,n_a);
%
figure;
%plot3(X,Y,Z,'r.','MarkerSize',5);
%hold on;
fill3(faults_X,faults_Y,faults_Z,S,'LineWidth',1);
axis equal
set(gca,'XLim',[-20,20],'YLim',[-15,15],'ZLim',[-10,0]);
colorbar;
colormap('jet');

%length width depth dip strike east north strike-slip dip-slip opening
fault_M = [l_n*ones(1,n_a);w_n*ones(1,n_a);-Z;-90*ones(1,n_a);s*ones(1,n_a);X;Y;...
    S;zeros(1,n_a);zeros(1,n_a)];

%% synthetic data
n_x = 100;
n_y = 60;
n_o = n_x*n_y;
obs_x = linspace(-25,25,n_x);
obs_y = linspace(-15,15,n_y);

[obs_X,obs_Y] = meshgrid(obs_x,obs_y);
obs_XYZ = [obs_X(:),obs_Y(:),zeros(n_o,1)]';

mu = 30e9;
nu = 0.25;
U3d = disloc3d(fault_M,obs_XYZ,mu,nu);

%% spatial correlated noise
C = zeros();
for ii = 1:n_o
    for jj = 1:ii
        obs_r = sqrt((obs_XYZ(1,jj)-obs_XYZ(1,ii))^2 + (obs_XYZ(2,jj)-obs_XYZ(2,ii))^2);
        C(ii,jj) = 1.2*exp(-obs_r/2);
        C(jj,ii) = C(ii,jj);
    end
end
L = chol(C,'lower');

%%
az = 350;
incid = 40; % or 20
los_vector = [-sind(incid)*cosd(az),sind(incid)*sind(az),cosd(incid)];

obs_u = los_vector*U3d;
obs_n = 0.05*L*randn(n_o,1);
obs_u = obs_u(:) + obs_n;
obs_U = reshape(obs_u,n_y,n_x);


mask_U = normalize(abs(obs_U),'range');
mask_C = mask_U<0.9;
mask_c = mask_C(:);

figure;
imagesc(obs_x,obs_y,obs_U,'AlphaData',mask_C);
hold on;
plot(X,Y,'k-','LineWidth',3);
colorbar;
colormap(redblue);
axis xy;
wysiwyg(gcf,15,10);


%% green function
obs_G = zeros(n_o,n_a);
for ii = 1:n_a
    m_ii =  [l_n;w_n;-Z(ii);-90;s;X(ii);Y(ii);1;0;0];
    obs_G(:,ii) = los_vector*disloc3d(m_ii,obs_XYZ,mu,nu);
end

%% laplacian
laplacian = compute_laplacian(X(:),Y(:),Z(:),4);

%%
% obs_G = obs_G(mask_c,:);
% obs_u = obs_u(mask_c);

save('obs_G','obs_G');
save('obs_u','obs_u');
save('laplacian','laplacian');

function z = gaussian(x,y,xc,yc,sigma)
z = exp(-((x-xc).^2 + (y-yc).^2)./(2*sigma));
end

