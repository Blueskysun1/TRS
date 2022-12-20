% This is a main function that will invoke Algorithm KeyGen(PP), Algorithm Sign(PP, T, μ, SK)
% and Algorithm Verify(PP, T, μ, σ) to complete the implementation of the proposed scheme. 
% Considering the stability, the three functions mentioned above that are invoked will be 
% repeatedly executed 1000 times respectively. Therefore, we recommend the reader to try
% to test this program with a smaller ring size so that you do not wait too long.

clc; clear;
% Setup Algorithm outputs the public parameters PP. 

[n, m, d, q, A, f] = Setup;


% Algorithm KeyGen(PP) outputs (SK,PK).

tic;
% Timing tool: Start.
% Considering the instability of a single execution of KeyGen(PP) Algorithm by a personal computer,
% we repeat the execution of KeyGen(PP) Algorithm 1000 times.
for p = 1:1000
    [SK, PK] = KeyGen(m, d, n, A, f, q);
end
fprintf('The computational overhead of KeyGen(PP) Algorithm running 1000 times is %f sec.\n',toc);
% toc;
% Timing tool: End.


% In this module, the signer enters the ring size and the message to be signed.

% Let the public key set of the ring be Lpk, and to facilitate the implementation of
% the experiment, we let the index π of the signer be 1. For the public keys of other 
% ring members we can get them by random sampling in some specified set, which is 
% due to the determination M-LWE problem.
N = input("Please enter the size of the ring:");
% Lpk = randi([-q,q], N-1, d*n);
Lpk = [PK;randi([-q,q], N-1, d*n)];
% Input Message μ (μ is not a legal variable in matlab, thus we use miu to
% represent μ).
prompt = 'Please enter the message you need to sign(e.g. Hello World！): ';
miu = input(prompt, 's');


% Sign(PP, T, miu, SK) outputs σ = (θ, z, C).

% Note that d0 and B are not part of the ring signature in the output of the Sign
% function below. d0 and B are public and only facilitate the invocation of the 
% Verify algorithm.
tic;
% Timing tool: Start.
for p=1:1000
    [d0, B, theta, z, C] = Sign(q, n, d, m, A, f, Lpk, miu, SK);
end
% Timing tool: End.
fprintf('The computational overhead of Sign Algorithm running 1000 times is %f sec.\n',toc);


% Verify(PP, T, μ, σ) outputs 0/1.

tic;
% Timing tool: Start.
for v =1:1000
    [result] = Verify(q, n, d, m, A, f, Lpk, miu, d0, B, theta, z, C);
end
% Print verification result.
if (result)
    fprintf ('Ring signature is valid! Output 1\n' );
else
    fprintf ('Ring signature is invalid! Output 0\n' );
end
fprintf('The computational overhead of Verify Algorithm running 1000 times is %f sec.\n',toc);
% toc;
% Timing tool: End.
fprintf("\n");



