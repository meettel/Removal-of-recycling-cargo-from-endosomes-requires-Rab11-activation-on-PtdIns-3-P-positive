writerObj = VideoWriter('pippo');
writerObj.FrameRate = 5;
open(writerObj);           %% crea un filmato con le immagini segmentate
% clear gfret
for kk = 1:size(gfret_2_sonda_vps34in1_trf_001,1)
     
    frame=kk;
     
    figure('position', [0, 0, 770, 480]), title(['Activation vs distance from nucleus in frame ' num2str(kk) ' file 1 001']), hold on, plot(space,gfret_2_sonda_vps34in1_trf_001(kk,:),'bx-'), axis([0 300 0.9 1.5])

    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    marg = 30;
    rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    frame = getframe(gca,rect);
    ax.Units = 'normalized';
    frame = getframe(gcf,[10 30 700 450]);
%     s = size(frame.cdata);
% fprintf('%d %d\n', s(2), s(1))
    writeVideo(writerObj,frame);
    close all
end
close(writerObj);