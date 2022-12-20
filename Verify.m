function [result] = Verify(q, n, d, m, A, f, Lpk, miu, d0, B, theta, z, C)
    N = size(Lpk,1);
    % Compute di = d0 + i·θ, for all i∈[N].
    % Note that all di are stored in an N × (n × d) dimensional matrix.

    for i = 1:N
        di(i,1:n*d)=mod(d0+i*theta,q);
    end

    % Compute t0 = -A·z + ΣPK·c
    % Comments on the implementation of this operation can be found in KeyGen.m.
    pkc = zeros(n,2*d-1);
    sum = zeros(n,2*d-1);
    Q = zeros(n,d-1);
    R = zeros(n,2*d-1);
    for j=1:n
        for i = 0:n
            sum(j,:) = sum(j,:) + conv(-A(j,1+d*i:d+d*i),z(1,1+d*i:d+d*i));
        end
        sum(j,:) = sum(j,:) - [zeros(1,d-1),z(1,1+d*(j+n):d+d*(j+n))];
        for k = 1:N
            pkc(j,:) = pkc(j,:) + conv(Lpk(k,1+d*(j-1):d+d*(j-1)),C(k,:));
        end
        sum(j,:) = sum(j,:) + pkc(j,:);

        [Q(j,:),R(j,:)] = deconv(sum(j,:),f);
        t0(j,:) = mod(R(j,d:2*d-1),q);
    end

    % Compute t1 = -B·z + Σdi·c
    %The process is similar to the calculation of t0.
    dic = zeros(n,2*d-1);
    sum = zeros(n,2*d-1);
    Q = zeros(n,d-1);
    R = zeros(n,2*d-1);
    for j=1:n
        for i = 0:n
            sum(j,:) = sum(j,:) + conv(-B(j,1+d*i:d+d*i),z(1,1+d*i:d+d*i));
        end
        sum(j,:) = sum(j,:) - [zeros(1,d-1),z(1,1+d*(j+n):d+d*(j+n))];
        for k = 1:N
            dic(j,:) = dic(j,:) + conv(di(k,1+d*(j-1):d+d*(j-1)),C(k,:));
        end
        sum(j,:) = sum(j,:) + dic(j,:);
        [Q(j,:),R(j,:)] = deconv(sum(j,:),f);
        t1(j,:) = mod(R(j,d:2*d-1),q);
    end


    % Compute c = Σci, for all i∈[N].
    % All ci are stored in matrix C.
    c_zero = zeros(1,1);
    for k=1:N
        c_zero = c_zero + C(k,:);
    end                         
    % Make c and all ci in the same group of additive mod 3.
    % They are all in the same challenge space Dc.
    c = mod(c_zero,3);
    for i=1:d
        if(c(1,i) == 2)
            c(1,i) = -1;
        end
    end

    % To invoke the hash function, similarly, the input to the function
    % needs to be pre-processed first.
    t0_h = reshape(t0,1,d*n);
    t1_h = reshape(t1,1,d*n);
    miu_h = [double(miu), zeros(1,d*n-length(double(miu)))];
    Verify_c = hash ([Lpk; miu_h; t0_h; t1_h; theta],'SHA-512');
    % The result Verify_c of the hash output (character vector) is transformed
    % into a 1*128 vector hash_c composed of its ASCII code, which needs to
    % be assigned to c after ASCII mod 3 first, and then change 2 to -1.
    hash_c = [];
    for i=1:strlength(Verify_c)
        hash_c(1,i) = mod(Verify_c(1,i),3);
        if(hash_c(1,i) == 2)
            hash_c(1,i) = -1;
        end
    end

    % Determine whether the infinite norm of z exceeds md^2-d and c = hash_c
    % (Σci = hash ([[Lpk;t0_h;t1_h;theta],miu_h])).
    for i=1:m*d
        if(z(1,i) > m*d^2-d && z(1,i) < -m*d^2-d)
            break
        end
    end

    for j=1:d
        if(c(1,j) ~= hash_c(1,j))
            break
        end
    end
    if (i == m*d && j == d)
        result = 1;
    else
        result = 0;
    end
end