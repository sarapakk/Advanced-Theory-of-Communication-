function R_IYmat = FTD(R_IY, k1, k2)
    R_IY = R_IY(end:-1:1);
    R_IYmat  = zeros(k1+1, k2);
    for i = 1:k1+1
        R_IYmat(i, :) = R_IY(k2:-1:1);
        R_IY = circshift(R_IY, [-1, 0]);
    end
end