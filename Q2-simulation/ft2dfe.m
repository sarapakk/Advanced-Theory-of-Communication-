function [R_Ymat, R_Y] = ft2dfe (x, k1, SNR)
    XX = conv(x, x);
    
    X = zeros(size(XX));
    X(1:length(x)) = x;
    X = circshift(X, [0, (length(XX)-length(x))/2]);
    
    R_Y = zeros(1, 2*k1+1);
    R_Y(1:length(X)) = XX + X*1/SNR;
    R_Y = circshift(R_Y, [0, (length(R_Y)-length(XX))/2]);
    R_Y = R_Y(end:-1:1);
    R_Y = circshift(R_Y, [0, -floor(length(R_Y)/2)]);
    
    R_Ymat = zeros(k1+1); 
    for i = 1:k1+1
        R_Ymat(i,:) = R_Y(1:k1+1);
        R_Y = circshift(R_Y, [0, 1]);
    end
end