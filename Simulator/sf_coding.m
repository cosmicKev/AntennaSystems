%*************** Space-Frequency Coding ***************
function [tx1, tx2] = sf_coding(simb)
    tx1 = zeros(1,768);
    tx1(1:2:end) = simb(1:2:end)/sqrt(2);
    tx1(2:2:end) = -(conj(simb(2:2:end)))/sqrt(2);

    tx2 = zeros(1,768);
    tx2(1:2:end) =        simb(2:2:end)/sqrt(2);
    tx2(2:2:end) =   conj(simb(1:2:end))/sqrt(2);
end