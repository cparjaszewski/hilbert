function value = innerfun_nlo_g( argument, constValue, fun)
noArgs = size(argument);

if (abs(argument-constValue)<1E-10)
    value = 0;
else
    value = ((fun.*argument)./((argument+constValue).*(argument-constValue)));
end;

for i=1:noArgs(2)  
  if (abs(argument(i)-constValue)<1E-10)
     value(i) = 0;
  end;
end;

