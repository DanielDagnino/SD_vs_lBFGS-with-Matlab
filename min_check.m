%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
pos = [1, 21, 1920, 1060];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos1b = [1, 1920+1, 1920, 1080];
pos2b = [1, 1920+221, 1920, 860];
pos3b = [1, 1920+521, 1920, 560];
pos1 = [1, 1, 1920, 1080];
pos2 = [1, 221, 1920, 860];
pos3 = [1, 521, 1920, 560];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test.
% func = x^2 + (10*y)^2;

xsol = [0, 0];
x0 = [10, 1];
niter = 80;
step0 = 1;
mem = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SD.
step = step0;
x = x0;

xsaved = zeros(2,niter);
xsaved(:,1) = x; 
error_SD = zeros(1,niter);
% error_SD(1) = (xsol-x)*(xsol-x)';
error_SD(1) = x(1)^2 + (10*x(2))^2;
for iter=2:niter
  
  dir = -[ 2*x(1), 2*(10*x(2))*10 ];
  
  x_new = x + step*dir/sqrt(dir*dir');
  
%   error_new = (xsol-x_new)*(xsol-x_new)';
  error_new = x_new(1)^2 + (10*x_new(2))^2;
  
  if error_SD(iter-1) > error_new
    error_SD(iter) = error_new;
    x = x_new;
  else
    error_SD(iter) = error_SD(iter-1);
    step = step/2;
  end
  
  xsaved(:,iter) = x;  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lBFGS.
step = step0;
x = x0;

xsaved2 = zeros(2,niter);
xsaved2(:,1) = x; 
error_lBGFS = zeros(1,niter);
% error_lBGFS(1) = (xsol-x)*(xsol-x)';
error_lBGFS(1) = x(1)^2 + (10*x(2))^2;

grad = zeros(2,mem);
mod = zeros(2,mem);
mem_saved = 0;
for iter=2:niter
  
  grad_aux = [ 2*x(1), 2*(10*x(2))*10 ];
  
  for k=2:mem
    grad(:,k-1) = grad(:,k);
    mod(:,k-1) = mod(:,k);
  end
  grad(:,mem) = grad_aux;
  mod(:,mem) = x;
  
  mem_saved = min( mem_saved+1, mem );
  dir = lBFGS( 1, mem_saved, grad, mod )';
%   dir = -grad_aux;
  
  x_new = x + step*dir/sqrt(dir*dir');
  
%   error_new = (xsol-x_new)*(xsol-x_new)';
  error_new = x_new(1)^2 + (10*x_new(2))^2;
  
  if error_lBGFS(iter-1) > error_new
    error_lBGFS(iter) = error_new;
    x = x_new;
  else
    error_lBGFS(iter) = error_lBGFS(iter-1);
    step = step/2;
    mem_saved = 0;
  end
  
  xsaved2(:,iter) = x;  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot.
h = figure(1); clf(h); set(h,'OuterPosition',pos1); set(gcf,'color','w');

for niteraux=2:niter

    %--------------------------------------------------------------%
    subplot(2,2,1);
    p = plot( 1:niteraux, log(error_SD(1:niteraux)/error_SD(1)), '-ob' );
    set(gca,'linewidth',2,'fontsize',20);
    xlim([1 niteraux]); set(p,'linewidth',2);
    title('Log Error: Steepest Descent'); xlabel('Iteration','fontsize',12);

    subplot(2,2,2);
    p = plot( 1:niteraux, log(error_lBGFS(1:niteraux)/error_lBGFS(1)), '-og' );
    set(gca,'linewidth',2,'fontsize',20);
    xlim([1 niteraux]); set(p,'linewidth',2);
    title('Log Error: l-BFGS'); xlabel('Iteration','fontsize',12);
    
    %--------------------------------------------------------------%
    subplot(2,2,3);
    xx = xsaved(1,(1:niteraux));
    yy = xsaved(2,(1:niteraux));
    p = plot( xx, yy, '-b' ); hold on;
    plot( 0, 0, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'MarkerSize',5 );
    set(gca,'linewidth',2,'fontsize',20);
    xlim([1 niteraux]); set(p,'linewidth',2);
    title('Trajectory: Steepest Descent'); xlabel('x','fontsize',16); ylabel('y','fontsize',16);
    xlim([-x0(1) x0(1)]); ylim([-x0(2) x0(2)]);
    
    x = linspace(-10,10);
    y = linspace(-1,1);
    [X,Y] = meshgrid(x,y);
    Z = X.^2 + (10*Y).^2;
    p2 = contour(X,Y,Z,10,'--k'); hold on;
    
    subplot(2,2,4);
    xx = xsaved2(1,(1:niteraux));
    yy = xsaved2(2,(1:niteraux));
    p = plot( xx, yy, '-b' ); hold on;
    
    x = linspace(-10,10);
    y = linspace(-1,1);
    [X,Y] = meshgrid(x,y);
    Z = X.^2 + (10*Y).^2;
    p2 = contour(X,Y,Z,10,'--k'); hold on;
    
    plot( 0, 0, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'MarkerSize',5 );
    set(gca,'linewidth',2,'fontsize',20);
    xlim([1 niteraux]); set(p,'linewidth',2);
    title('Trajectory: l-BFGS'); xlabel('x','fontsize',16); ylabel('y','fontsize',16);
    xlim([-x0(1) x0(1)]); ylim([-x0(2) x0(2)]);

end





