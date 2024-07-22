function data=QpskModulation(bits)
    bits=reshape(bits,2,[]);
    data=(1-2*bits(1,:)+1j*(1-2*bits(2,:)))/sqrt(2);
end