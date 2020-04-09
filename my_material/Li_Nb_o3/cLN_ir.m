function n_c_e = cLN_ir(f,T)
% "Temperature-dependent Sellmeier equation for the index of refraction,
% ne, in congruent lithium niobate"
% Dieter H. Jundt 1997
lambda = 3e8./f/1e-6;
k = (T - 273 - 24.5) / ( T - 273 + 570.82);
n_c_e = sqrt( 4.5820 + ...
    ( 0.09921 + 5.2716e-8 * k ) ./ ( lambda.^2 - ( 0.21090 - 4.9143e-8 * k ).^2 ) + ...
    2.2971e-8 * k - 0.021940 * lambda.^2);
end