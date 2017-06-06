function hplot = plotResults(results, hplot)
% Plot raw data from dual angle video-oculography
%
% 2013
% Hannah Payne
% Raymond Lab, Stanford University
% hpayne@stanford.edu

if ~exist('hplot','var')
    hplot  = [];
end

if all(isnan(results.cr1a));   hplot = []; return;    end

if isempty(hplot)
    
    % Determine ylims
    alldata1 = [results.pupil1(:,1) results.cr1a(:,1) results.cr1b(:,1)];
    alldata2 = [results.pupil2(:,1) results.cr2a(:,1) results.cr2b(:,1)];
    ylims1 = [min(alldata1(:)) max(alldata1(:))];
    ylims2 = [min(alldata2(:)) max(alldata2(:))];
    
    
    % Plot camera 2
    h(1)=subplot(2,1,1);
    hplot(1,1) = line(results.time2, results.pupil2(:,1),'Color','k'); hold on;
    hplot(2,1) = line(results.time2, results.cr2a(:,1),'Color','r','LineStyle','--');
    hplot(3,1) = line(results.time2, results.cr2b(:,1),'Color','r');
    hplot(4,1) = line(results.time2, results.sync2,'Color','m','LineStyle','-');
    
    % Plot camera 1
    h(2)= subplot(2,1,2);
    hplot(1,2) = line(results.time1, results.pupil1(:,1),'Color','k'); hold on
    hplot(2,2) = line(results.time1, results.cr1a(:,1),'Color','r');
    hplot(3,2) = line(results.time1, results.cr1b(:,1),'Color','r','LineStyle','--');
    hplot(4,2) = line(results.time1, results.sync1,'Color','m','LineStyle','-');
    
    subplot(2,1,1);
    set(gca,'XColor','w')
    ylabel('Cam2 (px)')
    legend(hplot(:,1),{'pupilx','crx1''','crx2','sync'})
    ylims2 = [nanmean(results.pupil2(:,1))*.5 nanmean(results.pupil2(:,1))*1.5];
    ylim(ylims2);
    
    subplot(2,1,2);
    linkaxes(h,'x')
    xlabel('Time (s)')
    ylabel('Cam1 (px)')
    legend(hplot(:,2),{'pupilx','crx1','crx2''','sync'})
    ylims1 = [nanmean(results.pupil1(:,1))*.5 nanmean(results.pupil1(:,1))*1.5];
    ylim(ylims1);
    
else
    
    % Plot camera 2
    set(hplot(1,1), 'YData', results.pupil2(:,1));
    set(hplot(2,1), 'YData', results.cr2a(:,1));
    set(hplot(3,1), 'YData', results.cr2b(:,1));
    set(hplot(4,1), 'YData', results.sync2);
    
    % Plot camera 1
    set(hplot(1,2), 'YData', results.pupil1(:,1));
    set(hplot(2,2), 'YData', results.cr1a(:,1));
    set(hplot(3,2), 'YData', results.cr1b(:,1));
    set(hplot(4,2), 'YData', results.sync1);
    
    xlim([0 results.time1(find(~isnan(results.pupil1(:,1)),1,'last'))]);
end


success = 1;

