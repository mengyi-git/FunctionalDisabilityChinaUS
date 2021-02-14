function plotMeanCI_clhls_resid(xAge, y1, y2, y1Static, y2Static, y1Trend, y2Trend)

figure;
hold on

ciplot(quantile(y1, 0.025, 2), quantile(y1, 0.975, 2), xAge, rgb('light pink'))
plot(xAge, mean(y1, 2), 'r-', 'LineWidth', 1.5)
plot(xAge, y1Static, 'r:', 'LineWidth', 2)
% plot(xAge, y1Trend, 'r--', 'LineWidth', 2) % virtually overlaps with mean(y1, 2)

ciplot(quantile(y2, 0.025, 2), quantile(y2, 0.975, 2), xAge, rgb('light blue'))
plot(xAge, mean(y2, 2), 'b-', 'LineWidth', 1.5)
plot(xAge, y2Static, 'b:', 'LineWidth', 2)
% plot(xAge, y2Trend, 'b--', 'LineWidth', 2) % virtually overlaps with mean(y2, 2)

hold off

xlim([xAge(1), xAge(end)])
ylim([0, 1])

xlabel('Age')
ylabel('Probability')

%%
% plot(xAge, quantile(y1, 0.975, 2), 'b--', 'LineWidth', 1)
% plot(xAge, quantile(y2, 0.975, 2), 'r--', 'LineWidth', 1)
% plot(xAge, quantile(y1, 0.025, 2), 'b--', 'LineWidth', 1)
% plot(xAge, quantile(y2, 0.025, 2), 'r--', 'LineWidth', 1)