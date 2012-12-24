function value = innerfun_chi_g( imre, ws, w1)
noArgs = size(imre);

value = 1:noArgs(2);

for i=1:noArgs(2)
  if ((ws(i)-w1)<1E-20)
     value(i) = 0;
  else
     value(i) = imre(i) / (ws(i)^2-w1^2);
  end;
end;

