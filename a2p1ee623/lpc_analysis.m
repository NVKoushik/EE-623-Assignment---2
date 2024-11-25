function [a, G, residual] = lpc_analysis(frame, order)
    a = lpc(frame, order);
    residual = filter(a, 1, frame);
    G = sqrt(sum(residual.^2) / length(residual));
end
