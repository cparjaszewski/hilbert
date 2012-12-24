function res = ZScan_fit(opts)  
    % Read data from opts
    lambda = opts.lambda;
    filename = opts.filename; % TODO: validate file

    % READ data from file:
    [col1,col2,~,~] = textread(filename,'%f %f %f %f');
    
    nop = length(col1); % number of measured points
    consts = Consts(lambda, nop);

    zp0 = 1:consts.nop;
    for zz = 1:consts.nop 
        zp0(zz) = consts.zr*zz/consts.nop-consts.zr/2;
    end

    figure
    subplot(2,2,1), scatter(zp0,col2,5,'b','filled');
    title('Experimental data for 3mm, transmittance vs z');

    right_average = 0;
    left_average = 0;
    for zz=1:round(consts.nop/2)
        right_average = right_average+col2(zz+round(consts.nop/2)-1);
        left_average = left_average+col2(zz);
    end
    if left_average>right_average,         if_df0_negative = 1;
    else                                   if_df0_negative = 0;
    end

    scol2 = 0;
    for zz = 1:consts.nop, scol2 = scol2 + (col2(zz) - sum(col2)/consts.nop)^2; end
    sd = sqrt(scol2/(consts.nop-1));

    if abs(col2(1)-sum(col2)/consts.nop)>sd, col2(1)=sum(col2)/consts.nop; end
    
    if abs(col2(consts.nop)-sum(col2)/consts.nop)>sd, col2(consts.nop)=sum(col2)/consts.nop; end
    
    for zz = 2:(consts.nop-1)
        if abs(col2(zz)-col2(zz-1))+abs(col2(zz)-col2(zz+1))>2*sd
            col2(zz)=(col2(zz-1)+col2(zz+1))/2;
        end
    end

    subplot(2,2,2), scatter(zp0,col2,5,'b','filled');
    title(['Experimental data for ', opts.dataType, ' without sticking points']);

    for ii = 10:1:(consts.nop-10)
        for zz = -9:1:10
            frag(zz+10) = col2(ii+zz);
            xx(zz+10) = zz+10;
        end
        R_factor = corrcoef(xx,frag);
        linearity(ii-9) = abs(R_factor(1,2));
    end

    subplot(2,2,3), plot(consts.zr/consts.nop*((10-50):1:((consts.nop-10)-50)),linearity,'r');
    title('Linearity of the plot');

    right_range = 0;
    left_range = 0;
    for ii = round(length(linearity)/2):1:length(linearity)
        right_range = right_range + 1;
        if linearity(ii)<0.5
            break
        end
    end
    for ii = round(length(linearity)/2):1:length(linearity)
        jj = round(length(linearity))-ii+1;
        left_range = left_range + 1;
        if linearity(jj)<0.5
            break
        end
    end

    horizontal_shift = round((right_range-left_range)/2);

    col22=col2;
    for zz=(abs(horizontal_shift)+1):1:(consts.nop-(abs(horizontal_shift)+1))
        col22(zz)=col2(zz+horizontal_shift);
    end
    col2 = col22;

    new_range = (right_range+left_range)/2;
    beg = round(consts.nop/2-new_range)-10;
    fin = round(consts.nop/2+new_range)+10;

    w0e = sqrt(2*new_range/1.7/consts.k); % approximate w0 from experiment
    
    col2s = sort(col2);
    col2max = 0;
    col2min = 0;
    for zz=1:10
        col2max = col2max+col2s(length(col2s)-zz+1)/10;
        col2min = col2min+col2s(zz)/10;
    end
    zero_level = (col2max+col2min)/2;

    col2=col2/zero_level; % normalization to 1

    col2s = sort(col2);
    col2max = 0;
    col2min = 0;
    for zz=1:10
        col2max = col2max+col2s(length(col2s)-zz+1)/10;
        col2min = col2min+col2s(zz)/10;
    end
    dfe = (0-1)^if_df0_negative*(col2max-col2min)/0.406; % approximate delta phi from experiment
    
    subplot(2,2,4), scatter(zp0,col2,5,'b','filled');
    title('Experimental data for 3mm silica, shifted to the center and normalized');    

    
    [zp, best] = mainLoop(opts, w0e, consts, dfe);
    
    figure
    scatter(zp,col2,5,'b','filled');
    hold on
    plot(zp,best.Tr,'r');
    title('Best fit for 3mm silica');
    
    % get the results
    res = Results;
    res.approx_w0 = 1000*w0e; 
    res.approx_delta_phi_0 = dfe;
    res.best_w0 = 1000*best.w0;
    res.delta_phi_0 = best.df0; 
end

function [zp, best] = mainLoop(opts, w0e, consts, dfe)
    dif = 0:100;    
    dif(1)=1000;
    index = 0;
    
    hWaitBar = waitbar(0,['Please wait - fitting for ',opts.dataType]);
    for w0_index=1:10
        w0 = w0e+2*(w0_index-5)*10^(0-3); % [mm] waist of the beam at the focus
        z0 = consts.k*w0^2/2; % [mm] Rayleigh length
        wa = w0*sqrt(1+consts.d0^2/z0^2); % beam radius at the aperture plane

        dro = 3*wa/consts.ia;
        drc = consts.ra/consts.ia;

        for df0_index=1:10
            df0 = dfe+(df0_index-5)/10; % real part of nonlinear phase shift
            index = index+1;
            % Beginning of simulation code
            zp = 1:consts.nop; Tr = 1:consts.nop;
            for zz = 1:consts.nop
                z = consts.zr*zz/consts.nop-consts.zr/2+0.0000001*z0; % to avoid dividing by zero
                zp(zz)=z; % z for plot

                w = w0*sqrt(1+z^2/z0^2);
                df = -df0/(1+z^2/z0^2);
                d = consts.d0-z;
                R = z+z0^2/z; 
                g = 1+d/R;

                Tr(zz) = 0;
                PTo = 0;
                PTc = 0;  
                wm0 = 0:consts.ee; Rm=0:consts.ee; tm=0:consts.ee; wm = 0:consts.ee; fm = 0:consts.ee;
                for m = 0:consts.ee
                    wm0(m+1) = sqrt(w^2/(2*m+1));
                    dm = consts.k*wm0(m+1)^2/2;
                    wm(m+1) = wm0(m+1)*sqrt(g^2 + d^2/dm^2);
                    Rm(m+1) = d/(1 - g/(g^2 + d^2/dm^2));
                    tm(m+1) = atan(d/dm/g);             
                    product = 1;
                    for n = 0:m
                        product = product*(1+1i*(2*n-1)*consts.T/(4*pi));
                    end    
                    fm(m+1) = product*(1i*df)^m/factorial(m);
                end   
                for rr = 0:consts.ia % integration over radius
                    Eo = 0;
                    Ec = 0;
                    ro = 3*wa*rr/consts.ia;
                    rc = consts.ra*rr/consts.ia;
                    for m = 0:consts.ee
                        Eo = Eo + fm(m+1)*wm0(m+1)/wm(m+1)*exp(ro^2*(0-1/wm(m+1)^2-1i*consts.k/(2*Rm(m+1)))+1i*tm(m+1));
                        Ec = Ec + fm(m+1)*wm0(m+1)/wm(m+1)*exp(rc^2*(0-1/wm(m+1)^2-1i*consts.k/(2*Rm(m+1)))+1i*tm(m+1));
                    end
                    PTo = PTo + ro*dro*(abs(Eo))^2;
                    PTc = PTc + rc*drc*(abs(Ec))^2;
                end
                Tr(zz) = PTc/PTo;

            end

            % End of simulation code

            % Beginning of comparison-performing (fitting) code
            zero_level2 = (max(Tr)+min(Tr))/2;
            Tr = Tr/zero_level2;

%             dif(index+1) = 0;
%             for zz=beg:1:fin
%                 z = zr*zz/nop-zr/2+0.001*z0; % to avoid dividing by zero
%                 if(index+1) = dif(index+1) + abs((col2(zz)-Tr(zz)));
%             end

            if dif(index+1)==min(dif) && df0~=0
                best.Tr = Tr;
                best.w0 = w0;
                best.df0 = df0;
            end

            waitbar(index/100, hWaitBar);
        end
    end
    close(hWaitBar);
   
end