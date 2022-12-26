function [n, m, d, q, A, f] = Setup
% Set the parameters of M-LWE.
n = 7;
m = 15;
d = 128;
q = 67108753;

% Randomly sample a matrix A=[A0||In], where A0 denotes A' in the article,
% since A' cannot be defined as a legal variable in Matlab.
% Pick A0, we store the 7×8-dimensional polynomial vector in a 7*8*128 matrix.
A0 = randi([-(q-1)/2,(q-1)/2], n, d*(m-n));
I1 = zeros(1,d);
I1(1,d) = 1;
I0 = zeros(1,d);
I11 = [I1;I0;I0;I0;I0;I0;I0];
I12 = [I0;I1;I0;I0;I0;I0;I0];
I13 = [I0;I0;I1;I0;I0;I0;I0];
I14 = [I0;I0;I0;I1;I0;I0;I0];
I15 = [I0;I0;I0;I0;I1;I0;I0];
I16 = [I0;I0;I0;I0;I0;I1;I0];
I17 = [I0;I0;I0;I0;I0;I0;I1];
A = [A0,I11,I12,I13,I14,I15,I16,I17];  %matrix A=[A0||In]

% Constructing an irreducible polynomial f = x^d+1.
f = zeros(1,d+1);
f(1,1) = 1;
f(1,129) = 1;
end