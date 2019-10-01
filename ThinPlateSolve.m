% This is a function to calculate the deflections and stresses of a thin
% plate all the edges are simply supported under uniform distributed load.
% usage ThinPlateSolve001(0.5,0.5,0.01,70e9,0.3,0.5e6,40,40)
% usage ThinPlateSolve001(0.5,0.2,0.01,70e9,0.3,0.5e6,40,10)
% usage ThinPlateSolve001(0.2,0.5,0.01,70e9,0.3,0.5e6,10,40)
% usage ThinPlateSolve001(0.5,0.5,0.01,70e9,0.3,0.5e6,70,70)
a=0;
b=0;
height=0;
E=0;
mu=0;
qxy=0;
NDivX=0;
NDivY=0;

function ThinPlateSolve(a,b,height,E,mu,qxy,NDivX,NDivY) 
h = a/NDivX;
m = NDivX + 1;
k = b/NDivY;
n = NDivY + 1;
D = E*height*height*height/(12*(1-mu*mu));
r1 = 1/(h*h);
r2 = 1/(k*k);
dim = n*m;
A = zeros(dim,dim);
q = zeros(dim,1);
u = zeros(dim,1);
w = zeros(dim,1);
x_fig = zeros(m,1);
y_fig = zeros(n,1);
z_fig = zeros(n,m); %% for plotting the indices are reverse in matlab

Mx  = zeros(n,m);
My  = zeros(n,m);
Mxy = zeros(n,m);
Qx  = zeros(n,m);
Qy  = zeros(n,m);
SigmaX = zeros(n,m);
SigmaY = zeros(n,m);
TauXY  = zeros(n,m);

% Solving for the Poisson's equation grad*grad(u)= -q/D
% for boundary nodes parallel to x axis
    % ### For x = 0 ###
i=1;
for j = 1:n
    index = (j-1)*m + i;
    A(index,index)=1;
    q(index,1)=0;
    y_fig(j) = k*(j-1);
end
% ### For x = a ###
i = m;
for j=1:n
    index = (j-1)*m + i;
    A(index,index)=1;
    q(index,1)=0;
end

% ### for boundary nodes parallel to y axis ###

% ### For y = 0 ###
j=1;
for i=1:m
    index = (j-1)*m + i;
    A(index,index)=1;
    q(index,1)=0;
    x_fig(i) = h*(i-1);
end

% ### For y = b ###
j = n;
for i = 1 :m
    index = (j-1)*m + i;
    A(index,index)=1;
    q(index,1)=0;
end

% ### for interior nodes  ###
for i = 2 :m -1
    for j = 2 :n-1
        index_i_j = (j-1)*m + i;
        index_i_p_1_j = (j-1)*m + i + 1;
        index_i_m_1_j = (j-1)*m + i -1;
        index_i_j_p_1 = (j-1+1)*m + i;
        index_i_j_m_1 = (j-1-1)*m + i;
        A(index_i_j,index_i_j) = -2*(r1+r2);
        A(index_i_j,index_i_p_1_j) = r1;
        A(index_i_j,index_i_m_1_j) = r1;
        A(index_i_j,index_i_j_p_1) = r2;
        A(index_i_j,index_i_j_m_1) = r2;
        q(index_i_j,1)=-qxy/D;
    end
end
u = A\q;

% Solving for the Poisson's equation grad*grad(w)= u
% ### for boundary nodes parallel to y axis ###
% ### For y = 0 ### w = 0 -> q = 0
% ### For x = a ### w = 0 -> q = 0
% ### for boundary nodes parallel to y axis ###
% ### For y = 0 ### w = 0 -> q = 0
% ### For y = b ### w = 0 -> q = 0

% Coefficients for A and q will be the same for w as in u

% ### for interior nodes  ###
% Coefficients for A will be the same for w as in u but q will change
for i = 2 :m -1
    for j = 2 :n-1
        index_i_j = (j-1)*m + i;
        q(index_i_j,1)= u(index_i_j,1);
    end
end
w = A\q;
for i = 1 :m
    for j = 1 :n
        index = (j-1)* m + i ;
        z_fig(j,i) = w(index);
    end
end

% ### Calculation of Moments and Stresses ###
for i = 2:m-1
    for j = 2:n-1
        index_i_j = (j-1)*m + i;
        index_i_p_1_j = (j-1)*m + i + 1;
        index_i_m_1_j = (j-1)*m + i -1;
        index_i_j_p_1 = (j-1+1)*m + i;
        index_i_j_m_1 = (j-1-1)*m + i;
        index_i_p_1_j_p_1 = (j-1+1)*m + i + 1;
        index_i_p_1_j_m_1 = (j-1-1)*m + i + 1;
        index_i_m_1_j_p_1 = (j-1+1)*m + i -1;
        index_i_m_1_j_m_1 = (j-1-1)*m + i -1;
        Mx(j,i) = -D * ((w( index_i_p_1_j)-2*w(index_i_j)+w(index_i_m_1_j))/(h*h)+ mu * (w(index_i_j_p_1)-2*w(index_i_j)+w(index_i_j_m_1))/(k*k));
        My(j,i) = -D * ( mu*(w( index_i_p_1_j)-2*w(index_i_j)+ w(index_i_m_1_j))/(h*h) + (w(index_i_j_p_1)-2*w(index_i_j)+w(index_i_j_m_1))/(k*k));
        Mxy(j,i) = -D*(1-mu)*(w(index_i_p_1_j_p_1)-w(index_i_p_1_j_m_1)-w(index_i_m_1_j_p_1)+w( index_i_m_1_j_m_1))/(4*h*k);
        SigmaX(j,i) = -6 * Mx(j,i)/(height*height);
        SigmaY(j,i)= -6 * My(j,i)/(height*height);
        TauXY(j,i)= -6 * Mxy(j,i)/(height*height);
    end
end

% ### Calculation of Shear Forces ###
for i = 2 : m-1
    for j = 1: n
        index_i_p_1_j = (j-1)*m + i + 1;
        index_i_m_1_j = (j-1)*m + i -1;
        Qx(j,i) = -D * ( u(index_i_p_1_j)-u(index_i_m_1_j))/(2*h);
    end
end

for i = 1 :m
    for j = 2:n-1
        index_i_j_p_1 = (j-1+1)*m + i;
        index_i_j_m_1 = (j-1-1)*m + i;
        Qy(j,i) = -D * ( u(index_i_j_p_1)-u(index_i_j_m_1))/(2*k);
    end
end

index_1_1 = (1-1)*m + 1;
index_m_1 = (1-1)*m + m;
index_3_1 = (1-1)*m + 3;
index_2_1 = (1-1)*m + 2;
index_2_n = (n-1)*m + 2;
index_3_n = (n-1)*m + 3;
index_4_1 = (1-1)*m + 4;
index_4_n = (n-1)*m + 4;
index_m_m_1_1 = (1-1)*m + m-1;
index_m_m_2_1 = (1-1)*m + m-2;
index_m_m_3_1 = (1-1)*m + m-3;
index_m_m_1_n = (n-1)*m + m-1;
index_m_m_2_n = (n-1)*m + m-2;
index_m_m_3_n = (n-1)*m + m-3;
index_1_3 = (3-1)*m + 1;
index_1_2 = (2-1)*m + 1;
index_m_2 = (2-1)*m + m;
index_m_3 = (3-1)*m + m;
index_1_4 = (4-1)*m + 1;
index_m_4 = (4-1)*m + m;
index_1_n = (n-1)*m + 1;
index_1_n_m_1 = (n-1-1)*m + 1;
index_1_n_m_2 = (n-1-2)*m + 1;
index_1_n_m_3 = (n-1-3)*m + 1;
index_m_n   = (n-1)*m + m;
index_m_n_m_1 = (n-1-1)*m + m;
index_m_n_m_2 = (n-1-2)*m + m;
index_m_n_m_3 = (n-1-3)*m + m;
Qx(1,1) = -D*(1/(k*k*h)*(-2*w(index_2_1))+1/(h*h*h)*(w(index_3_1)-2*w(index_2_1)));
Qx(n,1) = -D*(1/(k*k*h)*(-2*w(index_2_n))+1/(h*h*h)*(w(index_3_n)-2*w(index_2_n)));
Qx(1,2) = -D*(1/(k*k*h)*(w(index_1_1)-w(index_3_1))+1/(2*h*h*h)*(w(index_4_1)-2*w(index_3_1)+2*w(index_1_1)+w(index_2_1)));
Qx(n,2) = -D*(1/(k*k*h)*(w(index_1_n)-w(index_3_n))+1/(2*h*h*h)*(w(index_4_n)-2*w(index_3_n)+2*w(index_1_n)+w(index_2_n)));
Qx(1,m-1) =  -D*(1/(k*k*h)*(w(index_m_m_2_1)-w(index_m_1))+1/(2*h*h*h)*(-w(index_m_m_1_1)-2*w(index_m_1)+2*w(index_m_m_2_1)-w(index_m_m_3_1)));
Qx(n,m-1) =  -D*(1/(k*k*h)*(w(index_m_m_2_n)-w(index_m_n))+1/(2*h*h*h)*(-w(index_m_m_1_n)-2*w(index_m_n)+2*w(index_m_m_2_n)-w(index_m_m_3_n)));
Qx(1,m) = -D*(1/(k*k*h)*(2*w(index_m_m_1_1))+1/(h*h*h)*(2*w(index_m_m_1_1)-w(index_m_m_2_1)));
Qx(n,m) = -D*(1/(k*k*h)*(2*w(index_m_m_1_n))+1/(h*h*h)*(2*w(index_m_m_1_n)-w(index_m_m_2_n)));
Qy(1,1) = -D*(1/(h*h*k)*(-2*w(index_1_2))+1/(k*k*k)*(w(index_1_3)-2*w(index_1_2)));
Qy(1,m) = -D*(1/(h*h*k)*(-2*w(index_m_2))+1/(k*k*k)*(w(index_m_3)-2*w(index_m_2)));
Qy(2,1) = -D*(1/(h*h*k)*(w(index_1_1)-w(index_1_3))+1/(2*k*k*k)*(w(index_1_4)-2*w(index_1_3)+2*w(index_1_1)+w(index_1_2)));
Qy(2,m) = -D*(1/(h*h*k)*(w(index_m_1)-w(index_m_3))+1/(2*k*k*k)*(w(index_m_4)-2*w(index_m_3)+2*w(index_m_1)+w(index_m_2)));
Qy(n-1,1) =  -D*(1/(h*h*k)*(w(index_1_n_m_2)-w(index_1_n))+1/(2*k*k*k)*(-w(index_1_n_m_1)-2*w(index_1_n)+2*w(index_1_n_m_2)-w(index_1_n_m_3)));
Qy(n-1,m) =  -D*(1/(h*h*k)*(w(index_m_n_m_2)-w(index_m_n))+1/(2*k*k*k)*(-w(index_m_n_m_1)-2*w(index_m_n)+2*w(index_m_n_m_2)-w(index_m_n_m_3)));
Qy(n,1) = -D*(1/(h*h*k)*(2*w(index_1_n_m_1))+1/(k*k*k)*(2*w(index_1_n_m_1)-w(index_1_n_m_2)));
Qy(n,m) = -D*(1/(h*h*k)*(2*w(index_m_n_m_1))+1/(k*k*k)*(2*w(index_m_n_m_1)-w(index_m_n_m_2)));

for j = 2 :n-1
    index_3_j = (j-1)*m + 3;
    index_2_j = (j-1)*m + 2;
    index_2_j_p_1 = (j-1+1)*m + 2;
    index_2_j_m_1 = (j-1-1)*m + 2;
    index_m_m_2_j = (j-1)*m + m-2;
    index_m_m_1_j = (j-1)*m + m-1;
    index_m_m_1_j_p_1 = (j-1+1)*m + m-1;
    index_m_m_1_j_m_1 = (j-1-1)*m + m-1;
    Qx(j,1) = -D*(1/(h*h*h)*(w(index_3_j)-2*w(index_2_j))+1/(h*k*k)*(w(index_2_j_p_1)-2*w(index_2_j)+w(index_2_j_m_1)));
    Qx(j,m) = -D*(1/(h*h*h)*(2*w(index_m_m_1_j)-w(index_m_m_2_j))+1/(h*k*k)*(2*w(index_m_m_1_j)-w(index_m_m_1_j_p_1)-w(index_m_m_1_j_m_1)));
end

for i = 3 :m-2
    index_i_p_2_1 = (1-1)*m + i+2;
    index_i_m_2_1 = (1-1)*m + i-2;
    index_i_p_1_1 = (1-1)*m + i+1;
    index_i_m_1_1 = (1-1)*m + i-1;
    index_i_p_2_n = (n-1)*m + i+2;
    index_i_m_2_n = (n-1)*m + i-2;
    index_i_p_1_n = (n-1)*m + i+1;
    index_i_m_1_n = (n-1)*m + i-1;
    Qx(1,i) = -D*(1/(2*h*h*h)*(w(index_i_p_2_1)-2*w(index_i_p_1_1)+2*w(index_i_m_1_1)-w(index_i_m_2_1))+1/(h*k*k)*(w(index_i_m_1_1)-w(index_i_p_1_1)));
    Qx(n,i) = -D*(1/(2*h*h*h)*(w(index_i_p_2_n)-2*w(index_i_p_1_n)+2*w(index_i_m_1_n)-w(index_i_m_2_n))+1/(h*k*k)*(w(index_i_m_1_n)-w(index_i_p_1_n)));
end

for i = 2 :m-1
    index_i_3 = (3-1)*m + i;
    index_i_2 = (2-1)*m + i;
    index_i_p_1_2 = (2-1)*m + i+1;
    index_i_m_1_2 = (2-1)*m + i-1;
    index_i_n_m_2 = (n-1-2)*m + i;
    index_i_n_m_1 = (n-1-1)*m + i;
    index_i_p_1_n_m_1 = (n-1-1)*m + i+1;
    index_i_m_1_n_m_1 = (n-1-1)*m + i-1;
    Qy(1,i) = -D*(1/(k*k*k)*(w(index_i_3)-2*w(index_i_2))+1/(k*h*h)*(w(index_i_p_1_2)-2*w(index_i_2)+w(index_i_m_1_2)));
    Qy(n,i) = -D*(1/(k*k*k)*(2*w(index_i_n_m_1)-w(index_i_n_m_2))+1/(k*h*h)*(2*w(index_i_n_m_1)-w(index_i_p_1_n_m_1)-w(index_i_m_1_n_m_1)));
end

for j = 3 :n-2
    index_1_j_p_2 = (j-1+2)*m + 1;
    index_1_j_m_2 = (j-1-2)*m + 1;
    index_1_j_p_1 = (j-1+1)*m + 1;
    index_1_j_m_1 = (j-1-1)*m + 1;
    index_m_j_p_2 = (j-1+2)*m + m;
    index_m_j_m_2 = (j-1-2)*m + m;
    index_m_j_p_1 = (j-1+1)*m + m;
    index_m_j_m_1 = (j-1-1)*m + m;
    Qy(j,1) = -D*(1/(2*k*k*k)*(w(index_1_j_p_2)-2*w(index_1_j_p_1)+2*w(index_1_j_m_1)-w(index_1_j_m_2))+1/(k*h*h)*(w(index_1_j_m_1)-w(index_1_j_p_1)));
    Qy(j,m) = -D*(1/(2*k*k*k)*(w(index_m_j_p_2)-2*w(index_m_j_p_1)+2*w(index_m_j_m_1)-w(index_m_j_m_2))+1/(k*h*h)*(w(index_m_j_m_1)-w(index_m_j_p_1)));
end

f1 = figure;
surf(x_fig,y_fig,z_fig);
xlabel('x');
ylabel('y');
zlabel('w');

colormap(flipud(jet));

f2 = figure;
surf(x_fig,y_fig,Mx);
xlabel('x');
ylabel('y');
zlabel('Mx');

colormap(flipud(jet));

f3 = figure;
surf(x_fig,y_fig,My);
xlabel('x');
ylabel('y');
zlabel('My');

colormap(flipud(jet));

f4 = figure;
surf(x_fig,y_fig,Mxy);
xlabel('x');
ylabel('y');
zlabel('Mxy');

f5 = figure;
surf(x_fig,y_fig,SigmaX);
xlabel('x');
ylabel('y');
zlabel('SigmaX');

f6 = figure;
surf(x_fig,y_fig,SigmaY);
xlabel('x');
ylabel('y');
zlabel('SigmaY');

f7 = figure;
surf(x_fig,y_fig,TauXY);
xlabel('x');
ylabel('y');
zlabel('TauXY');

f8 = figure;
surf(x_fig,y_fig,Qx);
xlabel('x');
ylabel('y');
zlabel('Qx');

f9 = figure;
surf(x_fig,y_fig,Qy);
xlabel('x');
ylabel('y');
zlabel('Qy');
end

ThinPlateSolve(0.5,0.5,0.01,70e9,0.3,0.5e6,70,70);