function SiO2_ir(f)
lambda=1e6.*my_c.c./f; %wavelength in micrometer
n_square=1+(0.6961663*lambda.^2)./(lambda.^2-0.0684043^2)+(0.4079426.*lambda.^2)./(lambda.^2-0.1162414^2)+(0.8974794.*lambda.^2)./(lambda.^2-9.896161^2);
n=sqrt(n_square);
end