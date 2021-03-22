clear all
close all
clc


k0 = 2*pi;
M  = 400;
X0 = linspace(-2,2,M);
Z0 = zeros(1,M);
% x0 = 10;
% z0 = 10;
Xs_min = -1.5; % rectangle dimensions
Xs_max =  1.5;
Z_min  = -2.7;
Z_max  = -0.7;

x0 = 2; z0 = 0.7;
F=@(x,z) ((-1i/4).*besselh(0,2,k0.*sqrt((x-x0).^2+(z-z0).^2))).^2;

% -----------------------------------------------
% In order to create a rectangular mesh in Matlab
% we can make the following procedure
% -----------------------------------------------
lambda = 1;
numDiv1 = (Xs_max - Xs_min)/(lambda/5); %number of subdivisions
numDiv2 = (Z_max  - Z_min)/(lambda/5); %number of subdivisions
hx=(Xs_max-Xs_min)/numDiv1;
hy=(Z_max -Z_min)/numDiv2;
Num = numDiv1*numDiv2;
xi  = Xs_min:hx:Xs_max; %all points in the x-axis
eta = Z_min:hy:Z_max;   %all points in the y-axis
[X,Z] = meshgrid(xi,eta); % create a grid of points
Y = F(X,Z); % function values
[elem,vertex] = surf2patch(X,Z,Y); % the main variables
numElem= numDiv1*numDiv2; %total number of elements
numVert= size(vertex,1); % total number of vertices
N=2; %number of Gauss points = NxN
G = zeros(M,Num);

for j = 1 : M
    
    x0 = X0(j); z0 = 0.7;
    F=@(x,z) ((-1i/4).*besselh(0,2,k0.*sqrt((x-x0).^2+(z-z0).^2))).^2;
    q = integral2(F,-1.5,1.5,-3.2,-1.2);
    
    integralTot=0;
    for i=1:numElem %compute the integral on each element
        v1=[vertex(elem(i,1),1),vertex(elem(i,1),2)];
        v2=[vertex(elem(i,2),1),vertex(elem(i,2),2)];
        v3=[vertex(elem(i,3),1),vertex(elem(i,3),2)];
        v4=[vertex(elem(i,4),1),vertex(elem(i,4),2)];
        vertices=[v1;v2;v3;v4];
        elemInteg=integQuad(F,vertices,N);
        integralTot=integralTot+elemInteg;
        G(j,i)  = elemInteg;
    end
    actualIntegVal(j) = (q); %exact known value
    errorInt(j) = abs(actualIntegVal(j)-integralTot); %absolute error
end



vnorm = vecnorm(G);
G = G*diag(1./vnorm);
[U,S,V]  =svd(G);

n=1; xs = 0; convg = false;
while ~convg
    sigma(n) = S(n,n);
    Res(n)  = sum(sigma.^2)/sum(diag(S).^2);
    
    if Res < 0.999
        n = n+1;
    else
        convg = true;
    end
end
Pr=conj(V(:,1:n));
G1 = G*Pr;
G1 = G;
selected = [];
available = 1:size(G1,1); % N = size(B1,1)

%% Select the first n-1 sensing locations
phi = [];
P = eye(size(G1,2));
for k = 1:n-1
    ip = zeros(1,length(available));
    for i = 1:length(available)
        ip(i) = norm(P*G1(available(i),:)');
    end
    [~,maximum] = max(ip);
    phi = [phi;G1(available(maximum),:)];
    selected = [selected, available(maximum)];
    available(maximum)=[];
    R = orth(phi');
    P = eye(size(G1,2)) - R*R';
end

%% Determine the remaining sensing locations
lam = -10;
tol = 0.12;
while (lam<tol)
    
    [V,D] = eig(phi'*phi);
    % Determine multiplicity of minimum eigenvalue
    U = V(:,1);
    for i = 2:n
        if(D(i)~=D(1))
            break;
        end
        U = [U,V(i,:)];
    end
    P = U*U';
    ip = zeros(1,length(available));
    for i = 1:length(available)
        ip(i) = norm(P*G1(available(i),:)');
    end
    [~,maximum] = max(ip);
    phi = [phi;G1(available(maximum),:)];
    selected = [selected, available(maximum)];
    if length(selected) == 33
        break;
    end
    available(maximum)=[];
    lam = D(1);
end

locs(:,1) = X0(selected);
locs(:,2) = Z0(selected)+z0;

scatter(locs(:,1),locs(:,2),'fill');
hold on
rectangle('Position',[Xs_min Z_min (Xs_max-Xs_min) (Z_max-Z_min)],'linewidth',2)
grid on
xlabel('x'); ylabel('z');
xlim([Xs_min - 0.5 Xs_max + 0.5]);
zz0 = z0 +0.5;
ylim([Z_min-0.5 zz0])
set(gca,'fontsize',20)


% wavelength

addpath('/Users/samvad/Desktop/Frenel Data')
%% Forward
incfn  = @(rho) besselh(0,2,k0*rho);

%% Set parameters of simulation %
% wavelength
lambda   = 1;
k0        = 2*pi/lambda;
eps1      = 1.3;
k1        = sqrt(eps1)*k0;
gamma     = 0.5772; e = exp(1);


%% Object Initialization (square object)
dis = 20 ;                    % discretization lam/20;
Zx  = 1*lambda; Zy = 1*lambda; cx = 0; cy = -1.75*lambda; shape='rect';
[theta1, rho1, w1, N1]  = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

hold on
scatter(rho1.*cos(theta1),rho1.*sin(theta1));


%% Matrix Initialization
n = N1;
A = zeros(2*n,2*n);
b = zeros(2*n,1); b1 = zeros(2*n,1);b2 = zeros(2*n,1);

% The two integrals are:
%oint[g1(r,p) grad(phi).n - grad(g1(r,p).n phi]dr = -phi_inc(p)  (1)
%oint[g2(r,p) grad(phi).n - grad(g2(r,p).n phi]dr = 0            (2)
gdiag1    = @(k,i) -1j/4*(w1(i) - k^2*w1(i)^3/48 - 1j*(2*w1(i)/pi*(log(w1(i)*k/(4*e))+gamma)));

%% Dashboard
for src = 1 : size(locs,1)
    
    Rho = locs(src,:);
    for i=1:N1
        p = [rho1(i)*cos(theta1(i)) rho1(i)*sin(theta1(i))];
        if i~=N1
            p_tmp = [rho1(i+1)*cos(theta1(i+1)) rho1(i+1)*sin(theta1(i+1))];
        else
            p_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        p_ed = p_tmp - p;
        for j=1:N1
            r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))];
            if j~=N1
                r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
            else
                r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
            end
            r_ed = r_tmp - r;
            nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            if i==j
                A(i,j)       = gdiag1(k0,i);
                A(i,j+N1)    = 1/2.0;
                A(i+N1,j)    = gdiag1(k1,i);
                A(i+N1,j+N1) = -1/2.0;
            else
                A(i,j)       =  w1(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
                A(i,j+N1)    = -w1(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
                A(i+N1,j)    =  w1(j) * glquad(@(t)green(k1,p,0.5,p_ed,r,t,r_ed),2);
                A(i+N1,j+N1) = -w1(j) * glquad(@(t)gradgreen(k1,p,0.5,p_ed,r,t,r_ed,nhat1),2);
            end
        end
        
        b(i,src) =  incfn(norm(p+0.5*p_ed - Rho));
    end
    xcoe(:,src) = A\b(:,src);
end


Nm     = size(locs,1);
B_pred = zeros(Nm,2*N1);
for i = 1:Nm
    p = [locs(i,1) locs(i,2)];  % ith segment
    
    for j=1:N1 %segments of the contour integral
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))]; %start of jth segment
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        B_pred(i,j)    = B_pred(i,j)    - w1(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        B_pred(i,j+N1) = B_pred(i,j+N1) + w1(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        
    end

end
E_rx = diag(B_pred*xcoe);

G_in = G(selected,:);

contrast = G_in\E_rx;
c1 = reshape(contrast,[10,15]);
imagesc(abs(c1))
