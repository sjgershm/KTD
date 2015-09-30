function KTD_plot(mode)
    
    switch mode
        
        case 'overshadowing_SOC'
            subplot(1,2,1);
            for i=1:2; results(i) = KTD_sim(i); end
            W = [results(1).W results(2).W];
            W = W(32:end,:);
            bar(W(1,:)); colormap bone
            set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'OV-A' 'OV-B'},'YLim',[0 0.025]);
            ylabel('Value','FontSize',25);
            xlabel('Condition','FontSize',25);
            title('Kalman TD','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,2);
            param = KTD_defparam;
            param.TD = 1;
            for i=1:2; results(i) = KTD_sim(i,param); end
            W = [results(1).W results(2).W];
            W = W(32:end,:);
            bar(W(1,:)); colormap bone
            set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'OV-A' 'OV-B'},'YLim',[0 0.025]);
            ylabel('Value','FontSize',25);
            xlabel('Condition','FontSize',25);
            title('TD','FontSize',25,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 1000 400]);
            
            figure;
            subplot(2,2,1);
            c = results(1).C(1:10,23);
            plot(c,'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            set(gca,'FontSize',20,'XLim',[0 size(c,1)+1],'XTick',1:size(c,1),'YLim',[-1 0.1]);
            xlabel('Phase 1 trial','FontSize',25);
            ylabel('Posterior covariance (AX)','FontSize',25);
            
            subplot(2,2,2);
            c = [results(1).K(1:10,3) results(2).K(1:10,23)];
            plot(c(:,1),'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            set(gca,'FontSize',20,'XLim',[0 size(c,1)+1],'XTick',1:size(c,1),'YLim',[-0.1 0.1]);
            xlabel('Phase 1 trial','FontSize',25);
            ylabel('Kalman gain (X)','FontSize',25);
            
            subplot(2,2,3);
            c = results(1).C(21:30,23);
            plot(c,'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            set(gca,'FontSize',20,'XLim',[0 size(c,1)+1],'XTick',1:size(c,1),'YLim',[-1 0.1]);
            xlabel('Phase 2 trial','FontSize',25);
            ylabel('Posterior covariance (AX)','FontSize',25);
            
            subplot(2,2,4)
            c = [results(1).K(21:30,3) results(2).K(21:30,23)];
            plot(c(:,1),'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k'); hold on;
            plot(c(:,2),'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','w');
            set(gca,'FontSize',20,'XLim',[0 size(c,1)+1],'XTick',1:size(c,1),'YLim',[-0.1 0.1]);
            xlabel('Phase 2 trial','FontSize',25);
            ylabel('Kalman gain (X)','FontSize',25);
            legend({'OV-A' 'OV-B'},'FontSize',25);
            
        case 'serial_compound_extinction'
            TD = [0 1];
            L = {'Kalman TD' 'TD'};
            for i = 1:2
                param = KTD_defparam;
                param.TD = TD(i);
                results = KTD_sim(3,param);
                W = [results(1).W results(2).W];
                subplot(1,2,i);
                bar(W(end,:)); colormap bone;
                set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'Ext' 'No-Ext'},'YLim',[0 0.15]);
                ylabel('Value','FontSize',25);
                xlabel('Group','FontSize',25);
                title(L{i},'FontSize',25,'FontWeight','Bold')
            end
            set(gcf,'Position',[200 200 1000 400]);
            
        case 'serial_compound_LI'
            TD = [0 1];
            L = {'Kalman TD' 'TD'};
            for i = 1:2
                param = KTD_defparam;
                param.TD = TD(i);
                results = KTD_sim(4,param);
                W = [results(1).W results(2).W];
                subplot(1,2,i);
                bar(W(end,:)); colormap bone;
                set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'Pre' 'No-Pre'},'YLim',[0 0.15]);
                ylabel('Value','FontSize',25);
                xlabel('Group','FontSize',25);
                title(L{i},'FontSize',25,'FontWeight','Bold')
            end
            set(gcf,'Position',[200 200 1000 400]);
            
        case 'serial_compound_LI_variance'
            results = KTD_sim(4); results = results(1);
            trial_length = 10;
            for n=1:length(results.model); k(n,:)=results.model(n).K([3 11]); end
            K = [k(3:trial_length:end,1) k(1:trial_length:end,2)]; K = K(1:10,:);
            for n=1:length(results.model); c(n,1)=results.model(n).C(3,3); end
            C = c(1:trial_length:end); C = C(1:10);
            subplot(1,2,1);
            plot(C,'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            set(gca,'FontSize',20,'XLim',[0 length(C)+1],'XTick',1:length(C));
            xlabel('Pre-exposure trial','FontSize',25);
            ylabel('Posterior variance (X)','FontSize',25);
            subplot(1,2,2);
            plot(K(:,1),'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            set(gca,'FontSize',20,'XLim',[0 size(K,1)+1],'XTick',1:size(K,1),'YLim',[0 0.3]);
            xlabel('Pre-exposure trial','FontSize',25);
            ylabel('Kalman gain (X)','FontSize',25);
            set(gcf,'Position',[200 200 1000 400]);
            
        case 'serial_compound_extinction_variance'
            results = KTD_sim(3); results = results(1);
            trial_length = 10;
            for n=1:length(results.model); k(n,:)=results.model(n).K(13); end
            K = k(3:trial_length:end); K = K(11:end);
            for n=1:length(results.model); c(n,1)=results.model(n).C(3,11); end
            C = c(1:trial_length:end); C = C(1:10);
            %             subplot(1,2,1);
            plot(C,'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            set(gca,'FontSize',20,'XLim',[0 length(C)+1],'XTick',1:length(C),'YLim',[-0.01 0.1]);
            xlabel('Conditioning trial','FontSize',25);
            ylabel('Posterior covariance (ZX)','FontSize',25);
            %             subplot(1,2,2);
            %             plot(K,'-ok','LineWidth',4,'MarkerSize',10,'MarkerFaceColor','k');
            %             set(gca,'FontSize',20,'XLim',[0 size(K,1)+1],'XTick',1:size(K,1),'YLim',[-0.3 0.3]);
            %             xlabel('Extinction trial','FontSize',25);
            %             ylabel('Kalman gain (Z)','FontSize',25);
            %             set(gcf,'Position',[200 200 1000 400]);
            %
            %         case 'latent_inhibition'
            %             param = KTD_defparam;
            %             results = KTD_sim(4,param);
            %             subplot(1,2,1);
            %             bar(results.W(end,:)); colormap bone
            %             set(gca,'FontSize',20,'XLim',[0.5 4.5],'XTick',1:4,'XTickLabel',{'Pre-B' 'Pre-A' 'Pre-A+Sco' 'Pre-A+Sco2'},'YLim',[0 0.5]);
            %             ylabel('Value','FontSize',25);
            %             xlabel('Condition','FontSize',25);
            %             title('KTD','FontSize',25,'FontWeight','Bold');
            %
            %             param.TD = 1;
            %             results = KTD_sim(4,param);
            %             subplot(1,2,2);
            %             bar(results.W(end,:)); colormap bone
            %             set(gca,'FontSize',20,'XLim',[0.5 4.5],'XTick',1:4,'XTickLabel',{'Pre-B' 'Pre-A' 'Pre-A+Sco' 'Pre-A+Sco2'},'YLim',[0 0.5]);
            %             ylabel('Value','FontSize',25);
            %             xlabel('Condition','FontSize',25);
            %             title('TD','FontSize',25,'FontWeight','Bold');
            %
            %             set(gcf,'Position',[200 200 1000 400]);
            
        case 'LI'
            results = KTD_sim(6);
            n = 10;
            figure;
            subplot(1,2,1);
            plot(results.W(n+1:end,1),'-ko','LineWidth',4,'MarkerFaceColor','k','MarkerSize',10); hold on;
            plot(results.W(n+1:end,2),'-ko','LineWidth',4,'MarkerFaceColor','w','MarkerSize',10);
            set(gca,'FontSize',20,'XLim',[0 n+1],'XTick',1:n,'YLim',[-0.1 1]);
            xlabel('Conditioning trial','FontSize',25);
            ylabel('Reward expectation','FontSize',25);
            legend({'Pre' 'No-Pre'},'FontSize',25,'Location','SouthEast');
            mytitle('A','left','FontSize',30,'FontWeight','Bold');
            subplot(1,2,2);
            plot(results.k(1:n),'-ko','LineWidth',4,'MarkerFaceColor','k','MarkerSize',10); hold on;
            set(gca,'FontSize',20,'XLim',[0 n+1],'XTick',1:n);
            xlabel('Pre-exposure trial','FontSize',25);
            ylabel('Kalman gain','FontSize',25);
            mytitle('B','left','FontSize',30,'FontWeight','Bold');
            set(gcf,'Position',[200 200 1000 400]);
            
        case 'recovery'
            M = [7 8 9 10];
            L{1} = {'Overshadowed' 'Unovershadowed'};
            L{2} = {'Blocked' 'Unblocked'};
            L{3} = {'Overexpected' 'Unoverexpected'};
            L{4} = {'Inhibited' 'Uninhibited'};
            T = {'A' 'B' 'C' 'D'};
            for i = 1:length(M)
                results = KTD_sim(M(i));
                subplot(2,2,i);
                bar(results.W(:,2)); colormap bone
                if i < 4
                    set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',L{i},'YLim',[0 1]);
                else
                    set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',L{i},'YLim',[-1 0]);
                end
                ylabel('Reward expectation','FontSize',25);
                mytitle(T{i},'Left','FontSize',30,'FontWeight','Bold')
            end
            set(gcf,'Position',[200 200 1000 800]);
            %tightfig
            
        case 'serial_backward_blocking'
            results = KTD_sim(11);
            n = 10;
            plot(results(1).W(11:end),'-ko','LineWidth',4,'MarkerFaceColor','w','MarkerSize',10); hold on;
            plot(results(2).W(11:end),'-ko','LineWidth',4,'MarkerFaceColor','k','MarkerSize',10);
            legend({'A+; B?' 'B+; A?'},'FontSize',25)
            set(gca,'FontSize',20,'XLim',[0 n+1],'XTick',1:n);
            xlabel('Trial','FontSize',25);
            ylabel('Value','FontSize',25);
            
        case 'second_order_LI'
            TD = [0 1];
            L = {'Kalman TD' 'TD'};
            for i = 1:2
                param = KTD_defparam;
                param.TD = TD(i);
                results = KTD_sim(16,param);
                w = [results(1).W(end) results(2).W(end)]
                subplot(1,2,i);
                bar(w); colormap bone
                set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'Pre' 'No-Pre'},'YLim',[0 0.08]);
                title(L{i},'FontSize',25,'FontWeight','Bold');
            end
            set(gcf,'Position',[200 200 800 300]);
            
        case 'second_order_ext'
            TD = [0 1];
            L = {'Kalman TD' 'TD'};
            for i = 1:2
                param = KTD_defparam;
                param.TD = TD(i);
                results = KTD_sim(17,param);
                w = [results(1).W(end) results(2).W(end)];
                subplot(1,2,i);
                bar(w); colormap bone
                set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'Ext' 'No-Ext'},'YLim',[0 0.08]);
                title(L{i},'FontSize',25,'FontWeight','Bold');
                ylabel('Value','FontSize',25);
            end
            set(gcf,'Position',[200 200 800 300]);
            
        case 'overshadowing_recovery'
            TD = [0 1];
            L = {'Kalman TD' 'TD'};
            for i = 1:2
                param = KTD_defparam;
                param.TD = TD(i);
                results = KTD_sim(19,param);
                w = results(2).W(end,:);
                subplot(1,2,i);
                bar(w); colormap bone
                set(gca,'FontSize',20,'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'X' 'Y'},'YLim',[0 0.4]);
                ylabel('Value','FontSize',25);
                title(L{i},'FontSize',25,'FontWeight','Bold');
            end
            set(gcf,'Position',[200 200 800 300]);
            
    end