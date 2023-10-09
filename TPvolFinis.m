clear all
close all

%% Variables du problème

l_x = 1;
l_y = 1;

L = 1;
Rc = L/4;

n_x = 50;
n_y = 50;
n_t = 10000;

a = 1;
Um = 6; % amplitude du bouillon

pas_x = 2*l_x/(n_x-1); 
pas_y = 2*l_y/(n_y-1);
pas_t = 1/n_t;

a1 = a*pas_t/(pas_x*pas_x);
a2 = a*pas_t/(pas_y*pas_y);

a3 = -pas_t/pas_x;
a4 = -pas_t/pas_y;

xi=[-l_x:pas_x:l_x];
yi=[-l_y:pas_y:l_y];

%% Construction de la matrice A (diffusion seule)

% A = zeros(n_x*n_y);
% 
% for i = 1: n_x
% 
%     for j = 1: n_y
% 
%         if i > 1 && i < n_x && j > 1 && j < n_y     % cas du noeud intérieur
% 
%             A(bijection(i,j,n_x), bijection(i-1,j,n_x)) = a1;
%             A(bijection(i,j,n_x), bijection(i+1,j,n_x)) = a1;
% 
%             A(bijection(i,j,n_x), bijection(i,j,n_x)) = 1-2*(a1+a2);
% 
%             A(bijection(i,j,n_x), bijection(i,j-1,n_x)) = a2;
%             A(bijection(i,j,n_x), bijection(i,j+1,n_x)) = a2;
% 
%         end
%     end
% end

%% Construction de la matrice A (schéma linéaire-linéaire)

% A = zeros(n_x*n_y);
% 
% for i = 1: n_x
% 
%     for j = 1: n_y
% 
%         if i > 1 && i < n_x && j > 1 && j < n_y     % cas du noeud intérieur
% 
%             Ue = -Um*(j*pas_y-l_y)/(((i+1/2)*pas_x-l_x)^2+(j*pas_y-l_y)^2);
%             Uw = -Um*(j*pas_y-l_y)/(((i-1/2)*pas_x-l_x)^2+(j*pas_y-l_y)^2);
% 
%             Vn = Um*(i*pas_x-l_x)/((i*pas_x-l_x)^2+((j+1/2)*pas_y-l_y)^2);
%             Vs = Um*(i*pas_x-l_x)/((i*pas_x-l_x)^2+((j-1/2)*pas_y-l_y)^2);
% 
%             A(bijection(i,j,n_x), bijection(i-1,j,n_x)) = a1-a3*Uw/2;
%             A(bijection(i,j,n_x), bijection(i+1,j,n_x)) = a1+a3*Ue/2;
% 
%             A(bijection(i,j,n_x), bijection(i,j,n_x)) = 1-2*(a1+a2)+a3*(Ue-Uw)/2+a4*(Vn-Vs)/2;
% 
%             A(bijection(i,j,n_x), bijection(i,j-1,n_x)) = a2-a4*Vs/2;
%             A(bijection(i,j,n_x), bijection(i,j+1,n_x)) = a2+a4*Vn/2;
% 
%         end
%     end
% end

%% Construction de la matrice A (schéma linéaire-upwind)

A = zeros(n_x*n_y);

for i = 1: n_x

    for j = 1: n_y

        if i > 1 && i < n_x && j > 1 && j < n_y     % cas du noeud intérieur

            Ue = -Um*(j*pas_y-l_y)/(((i+1/2)*pas_x-l_x)^2+(j*pas_y-l_y)^2);
            Uw = -Um*(j*pas_y-l_y)/(((i-1/2)*pas_x-l_x)^2+(j*pas_y-l_y)^2);

            Vn = Um*(i*pas_x-l_x)/((i*pas_x-l_x)^2+((j+1/2)*pas_y-l_y)^2);
            Vs = Um*(i*pas_x-l_x)/((i*pas_x-l_x)^2+((j-1/2)*pas_y-l_y)^2);

            A(bijection(i,j,n_x), bijection(i-1,j,n_x)) = a1-a3*max(Uw,0); 
            A(bijection(i,j,n_x), bijection(i+1,j,n_x)) = a1+a3*min(Ue,0);

            A(bijection(i,j,n_x), bijection(i,j,n_x)) = 1-2*(a1+a2)+a3*(max(Ue,0)-min(Uw,0))+a4*(max(Vn,0)-min(Vs,0)); 

            A(bijection(i,j,n_x), bijection(i,j-1,n_x)) = a2-a4*max(Vs,0);  
            A(bijection(i,j,n_x), bijection(i,j+1,n_x)) = a2+a4*min(Vn,0);

        end
    end
end

%% Construction de la matrice A (schéma linéaire-quick)

% A = zeros(n_x*n_y);
% 
% for i = 1: n_x
% 
%     for j = 1: n_y
% 
%         if i > 2 && i < n_x-1 && j > 2 && j < n_y-1     % cas du noeud intérieur
% 
%             Ue = -Um*(j*pas_y-l_y)/(((i+1/2)*pas_x-l_x)^2+(j*pas_y-l_y)^2);
%             Uw = -Um*(j*pas_y-l_y)/(((i-1/2)*pas_x-l_x)^2+(j*pas_y-l_y)^2);
% 
%             Vn = Um*(i*pas_x-l_x)/((i*pas_x-l_x)^2+((j+1/2)*pas_y-l_y)^2);
%             Vs = Um*(i*pas_x-l_x)/((i*pas_x-l_x)^2+((j-1/2)*pas_y-l_y)^2);
% 
%             A(bijection(i,j,n_x), bijection(i-1,j,n_x)) = a1+a3*(max(Ue,0)*(-1/4+1/8)-max(Uw,0)*(-1/4+1/8)-min(Uw,0)*(-1-1/4)); 
%             A(bijection(i,j,n_x), bijection(i+1,j,n_x)) = a1+a3*(max(Ue,0)*(1/4+1/8)+min(Ue,0)*(1-4)-max(Uw,0)*(1/4+1/8));
%             A(bijection(i,j,n_x), bijection(i+2,j,n_x)) = A(bijection(i,j,n_x), bijection(i+2,j,n_x)) + a3*(min(Ue,0)*(-1/4+1/8));
%             A(bijection(i,j,n_x), bijection(i-2,j,n_x)) = A(bijection(i,j,n_x), bijection(i-2,j,n_x)) + a3*(min(Uw,0)*(-1/4-1/8));
% 
%             A(bijection(i,j,n_x), bijection(i,j,n_x)) = 1-2*(a1+a2)+a3*(1+max(Ue,0)*(-1/4)+min(Ue,0)*(-1/8-3/4)-1-max(Uw,0)*(-1/4)-min(Uw,0)*(3/4-1/8))+a4*(1+max(Vn,0)*(-1/4)+min(Vn,0)*(-1/8-3/4)-1-max(Vs,0)*(-1/4)-min(Vs,0)*(3/4-1/8)); 
% 
%             A(bijection(i,j,n_x), bijection(i,j-1,n_x)) = a2+a4*(max(Vn,0)*(-1/4+1/8)-max(Vs,0)*(-1/4+1/8)-min(Vs,0)*(-1-1/4));
%             A(bijection(i,j,n_x), bijection(i,j+1,n_x)) = a2+a4*(max(Vn,0)*(1/4+1/8)+min(Vn,0)*(1-4)-max(Vs,0)*(1/4+1/8));
%             A(bijection(i,j,n_x), bijection(i,j+2,n_x)) = A(bijection(i,j,n_x), bijection(i,j+2,n_x)) + a3*(min(Vn,0)*(-1/4+1/8));
%             A(bijection(i,j,n_x), bijection(i,j-2,n_x)) = A(bijection(i,j,n_x), bijection(i,j-2,n_x)) + a3*(min(Vs,0)*(-1/4-1/8));
% 
%         end
%     end
% end

%% Construction de teta

% CI sur teta (CI centrée)

% teta = zeros(n_x*n_y,n_t);
% 
% for i = 1: n_x
% 
%     for j = 1: n_y
% 
% 
% 
%         if (i*pas_x-l_x)*(i*pas_x-l_x)+(j*pas_y-l_y)*(j*pas_y-l_y) <= Rc*Rc   % -l_ car CI centrée
% 
%             teta(bijection(i,j,n_x),1) = 1;
% 
%         end
%     end
% end

% CI sur teta (CI non centrée)

teta = zeros(n_x*n_y,n_t);

for i = 1: n_x

    for j = 1: n_y



        if (i*pas_x-(l_x/2)-l_x)^2+(j*pas_y-l_y)^2 <= Rc*Rc   % -l_ car CI centrée

            teta(bijection(i,j,n_x),1) = 1;

        end
    end
end


for t = 1:n_t-1

    teta(:,t+1) = A*teta(:,t);

    eta = reshape(teta(:,t),n_x,n_y);
    surf(xi,yi,eta);
    axis([-l_x,l_x,-l_y,l_y,0,1]);
    drawnow
end









