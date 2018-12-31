function Rb = Rb_SP(ustar,D,P)

nu = 1.46e-5.*1013./P;

Sc = nu./D; %schmidt number
Rb = 5*Sc.^(2/3)./ustar;