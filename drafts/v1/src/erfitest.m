function Y = erfitest(Z)
    %% This function has been found online - the credits are described
    %% underneath. Maple has already built-in erf function for complex
    %% arguments
    
    % function y=erfi( imag(z) ) = imag(erf(z))
    % Returns imaginary part of the error function of imaginary z

    % ERFI.M 4/30/92 Hilmar Schlegel 
    % Internet: hshlgaii@mailszrz.zrz.tu-berlin.de (52.32/13.25)
    % Latest change: 5/6/92
%     if real(Z) && ~imag(Z), Y = erf(Z); 
%     elseif imag(Z) && ~real(Z), Y = erfimag(Z); 
%     else error (
%     end;
%     
% end
% 
% function Y = erfimag(Z)


     t=real(Z); % arg check
     
     if imag(Z),disp('*** Warning: Imaginary part not real in ERFI '),end;
     
     [nz,mz]=size(Z);
     t=t(:);
     s1=t;
     s2=t;
     Z=t.^2;
     new=(1:length(t))';
     n=1;
     
     while true
         t(new)=t(new).*Z(new)/n;
         s2(new)=s1(new)+t(new)/(2*n+1);
         new=find(s2~=s1);
         if isempty(new),break,end
         n=n+1;
         s1(new)=s2(new);
     end
     
     Y=reshape(2*s2/sqrt(pi),nz,mz);
end