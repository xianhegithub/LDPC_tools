function err_rate = ldpcDecode(RxSig, newH, rho_dB, Niter_sumprod, Ns, alg_decode, cw)

[~, num_c] = size(newH);
cw_hat = zeros(num_c, Ns);

if alg_decode == 0

    cw_hat = (sign(RxSig)+1)/2;

elseif alg_decode == 1

    for msg_idx = 1:Ns
        cw_hat(:, msg_idx) = sumProduct(RxSig(:, msg_idx), newH, rho_dB, Niter_sumprod);
    end

elseif alg_decode == 2

    for msg_idx = 1:Ns
        cw_hat(:, msg_idx) = sumProductLog(RxSig(:, msg_idx), newH, rho_dB, Niter_sumprod);
    end

elseif alg_decode == 3

    for msg_idx = 1:Ns
        cw_hat(:, msg_idx) = minSum(RxSig(:, msg_idx), newH, Niter_sumprod);
    end

elseif alg_decode == 4

    for msg_idx = 1:Ns
        cw_hat(:, msg_idx) = bitFlipping(RxSig(:, msg_idx), newH, Niter_sumprod);
    end

end
[~, err_rate] = bitErr(cw,cw_hat);