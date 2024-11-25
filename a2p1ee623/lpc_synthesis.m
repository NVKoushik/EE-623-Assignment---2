function synthesized = lpc_synthesis(residual, a, G)
    synthesized = filter(1, a, residual) * G;
end
