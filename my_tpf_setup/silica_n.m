function n=silica_n(f)
%I. H. Malitson. Interspecimen comparison of the refractive index of fused silica, J. Opt. Soc. Am. 55, 1205-1208 (1965)
lambda=1e6.*my_c.c./f;
a=[0.6961663, 0.4079426, 0.8974794];
b=[0.0684043,0.1162414,9.896161];
chi=0;
for iter_n=1:length(a)
chi=chi+a(iter_n).*lambda.^2./(lambda.^2-b(iter_n).^2);
end
n=sqrt(1+chi);

end