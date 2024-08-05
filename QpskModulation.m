function data=QpskModulation(bits)
    % modulates bits to QPSK complex amplitudes
    % `bits` must be even-length
    bits=reshape(bits,2,[]);
    data=(1-2*bits(1,:)+1j*(1-2*bits(2,:)))/sqrt(2);
end