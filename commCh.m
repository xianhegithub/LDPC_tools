function RxSig = commCh(TxSig, channel, rho_dB)

rho_lin = 10^(rho_dB/10);
if channel == 'AWGN'
    N0 = 1/rho_lin;
    GN = sqrt(N0)*randn(size(TxSig));
    RxSig = TxSig + GN;
end
