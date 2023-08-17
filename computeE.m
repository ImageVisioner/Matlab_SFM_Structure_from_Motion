function E = computeE(F,K1,K2)
F1 = F/F(3,3);
E = K2'*F1*K1;
E = E/E(3,3);

end