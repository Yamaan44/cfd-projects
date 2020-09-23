function y = ADRtridiag(main,lower,upper,f)
% Mark Holmes (2020), edit by Mary McGuinn & Ariana Wetzel 
% [ a(1) c(1)                       ] [ y(1) ]  [ f(1) ]
% [ b(2) a(2) c(2)                  ] [ y(2) ]  [ f(2) ]
% [      b(3) a(3) c(3)             ] [      ]  [      ]
% [      ...  ...  ...              ] [ ...  ]= [ ...  ]
% [          ...  ...  ...          ] [      ]  [      ]
% [            b(n-1) a(n-1) c(n-1) ] [y(n-1)]  [f(n-1)]
% [                    b(n)   a(n)  ] [ y(n) ]  [f(n)  ]
%
% f must be a vector (row or column) of length n
I = length(f);               % f vector, length n
v = zeros(I,1);  
y = v;
w = main(1);
y(1) = f(1)/w;
for i=2:I
  v(i-1) = upper(i-1)/w;
  w = main(i) - lower(i)*v(i-1);
  y(i) = ( f(i) - lower(i)*y(i-1) )/w;
end
for j=I-1:-1:1
  y(j) = y(j) - v(j)*y(j+1);
end

end