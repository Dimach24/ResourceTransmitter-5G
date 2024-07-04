function out_seq = dmDescrambling(in_seq, N_Cell_ID, L_max, block_index)
    %dmDescrambling Procedure of revererse scrambling
    % after demodulation [7.3.3.1, TS 38.211] 
    arguments
        in_seq (1,:) % input sequence (boolean matrix)
        N_Cell_ID (1,1)
        L_max (1,1) % maximum number of candidate SS/PBCH blocks in half frame [5, TS 38.213]
        block_index (1,1) % candidate SS/PBCH block index
    end
    
    %init
    A = length(in_seq);
    s = zeros(1,A);
    M = A;

    %determinaton of nu
    block_index = fliplr(dec2bin(block_index));
    if L_max == 4
        nu = [block_index(2) block_index(1)];
    else
        nu = [block_index(3) block_index(2) block_index(1)];
    end
    nu = bin2dec(num2str(nu));

    %determination of c
    x1 = zeros(1,2000);
    x2 = zeros(1,2000);
    x1(1) = 1;
    x1(2:31) = 0;
    x2(1:31) = fliplr(dec2bin(N_Cell_ID,31)); %c_init = N_Cell_ID
    for n = 1:2000
        x1(n+31) = mod(x1(n+3)+x1(n),2);
        x2(n+31) = mod(x2(n+3)+x2(n+2)+x2(n+1)+x2(n),2);
        n1 = 1:160;
        c(n1) = mod(x1(n1+1600)+x2(n1+1600),2);
    end

    %determination of s
    i = 0;
    j = 0;
    while i < A
        s(1+i) = c(1+j+nu*M);
        j = j+1;
        i = i+1;
    end
    
    %descrambling procedure
    out_seq = zeros (1,A);
    for i = 1:A
    out_seq(i) = mod(in_seq(i)+ s(i),2);
    end
end