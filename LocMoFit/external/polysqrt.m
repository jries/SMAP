function pol = polysqrt(p)
%POLYSQRT Find a square root of a polynomial.
%   POL = POLYVAL(P) returns a vector POL, if it exists, such that
%   conv(POL,POL) = P. P is a vector whose elements are the coefficients 
%   of a polynomial in descending powers.

    deg_p = length(p) - 1;
    if( mod(deg_p,2) )
        pol = [];
        return
    end
    
    n_var = deg_p/2 + 1;
    pol = zeros(1,n_var);
    
    pol(1) = sqrt(p(1));
   
    for pos=2:n_var
        indep_term = p(pos) - calc_delta(pol(2:pos-1));
        pol(pos) = indep_term/(2*pol(1));
    end
    
    test_pol = conv(pol,pol);
    res = max(abs(p-test_pol));
    if(res > 1e-6)
        pol = [];
    end
end
    
function res = calc_delta(t)   
	k = length(t);
	res = 0;
	for i=1:k
        res = res + t(i)*t(k-i+1);
	end
end