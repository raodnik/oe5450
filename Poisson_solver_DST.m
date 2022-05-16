clc;
clear;

% Domain dimensions and grid size

L = 10;
N = 500; 
f = zeros(N-1,N-1);
d = L/(N-1);
                     

% Sources locations

f((N)/2+1, ((N))/5 + 1) =  -1.00;
f((N)/2+1, (4*(N))/5 + 1) = 1.00;

% Poisson Equation solution using FFT

F = dst2(f);
U = zeros(N-1,N-1);

for i = 1:N-1
    for j = 1:N-1
        U(i,j) = (0.25*F(i,j))/(1 - (cos((pi*(i+j))/(2*(N))))*(cos((pi*(i-j))/(2*(N)))));
    end
end

s = idst2(U);
u = zeros(N+1,N+1);
u(2:N,2:N) = s;

figure
imagesc(u);
colorbar;
axis off;


figure
contour(u);

%End of the Program
