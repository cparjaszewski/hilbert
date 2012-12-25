function plotNKSpectra(tab)
    clf;
    plot(tab(1,:), tab(2,:), tab(1,:), tab(3,:));
    legend('n','k');
end