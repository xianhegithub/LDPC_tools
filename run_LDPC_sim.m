Ns = 1000;                 % number of symbols in channel simulation
rho_dB = 5;       % normalised SNR in dB
rho_dB_low = -2;
rho_dB_up = 10;
num_r = 50;               % number of rows in the parity check matrix H  
num_c = 2*num_r;                % number of columns for the 
rate=.5;                % ldpc code rate, used to calculate number of columns in H, N = M/rate
method = 0;                   % parity check matrix generation method: 
% 0 evencol : For each column, place 1s uniformly at random
% 1 evenboth: For each column and row, place 1s uniformly at random
no4cycle = 1;           % if 1, remove length-4 cycles
checksPerCol = 3;       % Number of 1s per column for LDPC matrix
enc_strategy = 2;       % encoding strategy {0 = First; 1 = Mincol; 2 = Minprod}
modForm = 'BPSK';
channel = 'AWGN';
Niter_sumprod = 5;
alg_decode = 0;         % 0 for 'uncoded', 1 for 'sumProduct', 2 for 'sumProductLog'
% 3 for 'sumProductLogSimple', 4 for 'sumProductHard'

fprintf('Generating parity check matrix H of size %d X %d, and %d 1s per column\n'...
    , num_r, num_c, checksPerCol);

H = ldpcHMatrix(num_r, num_c, method, no4cycle, checksPerCol);

% Hack to handle low-rank H.
while rank(H) < num_r  
    fprintf('Dealing with low rank H\n');
    H = ldpcMatrix(num_r, num_c, method, no4cycle, checksPerCol);
end

% generate one random message/symbol
msg = round(rand(num_r, Ns));

% different encoding matrix for each message

[checkBits, newH] = ldpcEncode(msg, H, enc_strategy);
cw = [checkBits; msg];

% loop over different SNRs
err_rate = cell(1,5);
for ii = 1:5
    err_rate{ii} = zeros(length(rho_dB_low:rho_dB_up));
end

for alg_decode = 0:4
    if alg_decode == 0
        fprintf('decoding uncoded\n');
    elseif alg_decode == 1
        fprintf('decoding using sum product\n');
    elseif alg_decode == 2
        fprintf('decoding using sum product log\n');
    elseif alg_decode == 3
        fprintf('decoding using the min-sum decoder\n');
    elseif alg_decode == 4
        fprintf('decoding using bit flipping \n');
    else
        fprintf('Unknown decoding method\n');
    end

    for rho_dB_loop = rho_dB_low:rho_dB_up

        rho_dB = rho_dB_loop;
        % try transmission in AWGN
        TxSig = codedMod(cw, modForm);
        RxSig = commCh(TxSig, channel, rho_dB);

        err_rate{alg_decode+1}(rho_dB + 3) = ldpcDecode(RxSig, newH, rho_dB, Niter_sumprod, Ns, alg_decode, cw);
        
        fprintf('BER %f at SNR %.1f dB\n', err_rate{alg_decode+1}(rho_dB + 3), rho_dB);
        
    end  %for rho_dB_loop
    
    semilogy(rho_dB_low:rho_dB_up, err_rate{alg_decode+1})
    hold on
    
end  % for alg_decode = 0:4

legend('uncoded', 'sumProd', 'sumProdLog', 'min-sum', 'bit flipping');
