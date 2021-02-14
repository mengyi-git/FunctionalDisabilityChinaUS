function plotMeanCI_clhls_hrs(xAge, y1, y2, y1Static, y2Static)

figure;
hold on

ciplot(quantile(y1, 0.025, 2), quantile(y1, 0.975, 2), xAge, rgb('light blue'))
plot(xAge, mean(y1, 2), 'b-', 'LineWidth', 1.5)
plot(xAge, y1Static, 'b:', 'LineWidth', 2)

ciplot(quantile(y2, 0.025, 2), quantile(y2, 0.975, 2), xAge, rgb('light pink'))
plot(xAge, mean(y2, 2), 'r-', 'LineWidth', 1.5)
plot(xAge, y2Static, 'r:', 'LineWidth', 2)

hold off

xlim([xAge(1), xAge(end)])
ylim([0, 1])

xlabel('Age')
ylabel('Probability')
