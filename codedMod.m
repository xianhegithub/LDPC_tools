function TxSig = codedMod(cw, modForm)

if modForm == 'BPSK'
    
    % BPSK modulation
    bpskMod = 2*cw - 1;
end

TxSig = bpskMod;