function plotData_2(summaryTable)
    % Plot concentration vs. time
    figure
    
    subplot(1,2,1), 
    hold on
    plot(summaryTable(:,1), summaryTable(:,2:4));
    grid
    legend('[S1]','[S2]','[S3]')
    title('Substrates')
    xlabel('time (s)')
    ylabel('concentration (nM)')
    hold off
    
    subplot(1,2,2),
    hold on
    plot(summaryTable(:,1), summaryTable(:,5:7));
    grid
    legend('[P1]','[P2]','[P3]')
    title('Products')
    xlabel('time (s)')
    ylabel('concentration (nM)')
    hold off
    
end