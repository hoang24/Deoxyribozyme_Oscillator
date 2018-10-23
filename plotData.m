function plotData(summaryTable)
    % Display data
    figure;
    hold on;
    grid;
    plot(summaryTable(:,1), summaryTable(:,2:end))
    title('Deoxyribozyme Oscillator');
    xlabel('time (s)');
    ylabel('concentation (nM)');
    legend('[S1]', '[S2]', '[S3]', '[P1]', '[P2]', '[P3]');
    hold off;
end