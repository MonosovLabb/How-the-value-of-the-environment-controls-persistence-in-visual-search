function textPrintMeans(cellObj, net, sing, noRewXs, rewXs, vertOffset, nrewColor, rewColor)
    currY = ylim;
    newHigh = currY(2);
    dist = currY(2) - currY(1);
    
    % Net printing
    % NR: 1/10 way down
    temp = currY(2) - dist/10;
    text(noRewXs(1), temp, ['p=' mat2str(net.p(1),3) '   r=' mat2str(net.rho(1),3)], 'Color', nrewColor);
    
    % R: 2/10 way down
    temp = currY(2) - 2*dist/10;
    text(noRewXs(1), temp, ['p=' mat2str(net.p(2),3) '   r=' mat2str(net.rho(2),3)], 'Color', rewColor);
    
%     % Sing printing
%     % 1/10 above each point
%     for i = 1:4
%         % nrew
%         temp = mean(cellObj{i,1});
%         temp = temp + dist/30 + vertOffset;
%         newHigh = max(newHigh, temp);
%         text(noRewXs(i), temp, ['n=' mat2str(sing.n(i,1))]);
%         
%         % rew
%         temp = mean(cellObj{i+1,2});
%         temp = temp + dist/30 + vertOffset;
%         newHigh = max(newHigh, temp);
%         text(rewXs(i), temp, ['n=' mat2str(sing.n(i+1,2))]);
%     end
    
    ylim([currY(1) newHigh]);
end