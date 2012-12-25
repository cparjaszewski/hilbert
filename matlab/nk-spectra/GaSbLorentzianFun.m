function val = GaSbLorentzianFun(x)
    amplitude = [323.7 650.8 66.28 226 253.1 128.6 299.9 244.9];
    center = [309.3 515.2 855 786.3 222.3 719.8 413.2 626.1];
    width = [69.35 67.12 21.08 49.25 53.95 48.74 70.68 62.33];

    val = 0;
    for n=1:8, val = val + amplitude(n).*width(n)./((x-center(n)).^2+width(n).^2); end
    val = (1/pi) .* val;

end