[na,~]=size(PRn);
[nb,~]=size(PTn);

figure()
for ii = 1:na
%     for jj = 1:na
        plot(f,squeeze(real(PRn(ii,1,:))),'LineWidth',1)
        hold on
%     end
end

figure()
for ii = 1:nb
%     for jj = 1:nb
        plot(f,squeeze(real(PTn(ii,1,:))),'LineWidth',1)
        hold on
%     end
end
% plot(f,squeeze(real(PTn(2,2,:))),'LineWidth',1)
total = abs(PRn(1,1,:))+abs(PRn(2,1,:))+abs(PRn(3,1,:))+abs(PTn(1,1,:))+abs(PTn(2,1,:))+abs(PTn(3,1,:));
plot(f,squeeze(total))


%% 
[na,~]=size(PR1);
[nb,~]=size(PT1);


figure()
for ii = 1:na
%     for jj = 1:na
        plot(f,squeeze(real(PR1(ii,:))),'LineWidth',1)
        hold on
%     end
end

figure()
for ii = 1:nb
%     for jj = 1:nb
        plot(f,squeeze(real(PT1(ii,:))),'LineWidth',1)
        hold on
%     end
end
% plot(f,squeeze(real(PTn(2,2,:))),'LineWidth',1)
% total = abs(PR1(35+1,:))+abs(PR1(35+2,:))+abs(PR1(35+3,:))+abs(PT1(1,:))+abs(PT1(2,:))+abs(PT1(3,:));
% total = abs(PR1(35+1,:))+abs(PR1(35+2,:))+abs(PT1(1,:))+abs(PT1(2,:));
% total = abs(PR1(35+1,:))+abs(PT1(1,:));

total = +(PR1(35+1,:))+(PR1(35+2,:))+(PR1(35+3,:))+(PT1(1,:))+(PT1(2,:))+(PT1(3,:));
total = (PR1(35+1,:))+(PR1(35+2,:))+(PT1(1,:))+(PT1(2,:));
total = (PR1(35+1,:))+(PT1(1,:));

% total = abs((PR1(1,:)))+abs((PT1(1,:)));
plot(f,squeeze(total))

%% 
[na,~]=size(PR2);
[nb,~]=size(PT2);


figure()
for ii = 1:round(na/2)
%     for jj = 1:na
        plot(f,squeeze(real(PR2(35+ii,:,:))),'LineWidth',1)
        hold on
%     end
end

figure()
for ii = 1:nb
%     for jj = 1:nb
        plot(f,squeeze(real(PT2(ii,:,:))),'LineWidth',1)
        hold on
%     end
end
% plot(f,squeeze(real(PTn(2,2,:))),'LineWidth',1)
total = abs(PR2(36,1,:))+abs(PR2(37,1,:))+abs(PR2(38,1,:))+abs(PT2(1,1,:))+abs(PT2(2,1,:))+abs(PT2(3,1,:));
plot(f,squeeze(total))