%%

clear all;
close all;
clc;
warning off;

%% Setting Parameters

No_Anchors = 4;                 % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
NetworkSize1 = 100;             % We consider a 100by100 area that the mobile can wander

Anchor_Node_Loc   = [0    0;    % set the anchor at 4 vertices of the region
    NetworkSize1           0;
    0           NetworkSize1;
    NetworkSize1 NetworkSize1];

Mobile_Node_Loc  = NetworkSize1*rand(No_MobileNode,2);  % Building a random location for the mobile node
Distance = zeros(No_Anchors,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : No_Anchors
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the Scenario

figure('Name','Mobile Node:100, Anchor Node:4, N/W Size:100','NumberTitle','off');
clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements

Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors,No_MobileNode)-1/2);
No_Iteration = 5;

% Initial Guess (Random Locatio)

MobileLoc_Est = NetworkSize1*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)-Anchor_Node_Loc(:,1))./Dist_Esti ...          % x-coordinate
            (MobileLoc_Est(m,2)-Anchor_Node_Loc(:,2))./Dist_Esti];                          % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',8,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location',...
    'Location','NorthEast')

% Compute the Root Mean Squred Error

Err1 = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/100;
axis([-0.1 1.1 -0.1 1.1]*NetworkSize1)

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, Anchor Node:4, N/W Size:100','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%% Setting Parameters

No_Anchors = 4;                 % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
NetworkSize2 = 200;             % We consider a 200by200 area that the mobile can wander

Anchor_Node_Loc   = [0     0;   % set the anchor at 4 vertices of the region
    NetworkSize2           0;
    0           NetworkSize2;
    NetworkSize2 NetworkSize2];

Mobile_Node_Loc  = NetworkSize2*rand(No_MobileNode,2);  % Building a random location for the mobile node

Distance = zeros(No_Anchors,No_MobileNode);
for m = 1 : No_MobileNode
    for n = 1 : No_Anchors
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the Scenario

figure('Name','Mobile Node:100, Anchor Node:4, N/W Size:200','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements

Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors,No_MobileNode)-1/2);
No_Iteration = 5;

% Initial Guess (Random Locatio)

MobileLoc_Est = NetworkSize2*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)-Anchor_Node_Loc(:,1))./Dist_Esti ...      % x-coordinate
            (MobileLoc_Est(m,2)-Anchor_Node_Loc(:,2))./Dist_Esti];                      % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',8,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');

legend('Anchor locations','Mobile true location','Mobile estimated location',...
    'Location','NorthEast')

% Compute the Root Mean Squred Error
Err2 = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/100;
axis([-0.1 1.1 -0.1 1.1]*NetworkSize2)

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, Anchor Node:4, N/W Size:200','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%% Setting Parameters

No_Anchors = 4;                 % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
NetworkSize3 = 300;             % We consider a 300by300 area that the mobile can wander

Anchor_Node_Loc   = [0    0;    % Set the anchor at 4 vertices of the region
    NetworkSize3           0;
    0           NetworkSize3;
    NetworkSize3 NetworkSize3];

Mobile_Node_Loc  = NetworkSize3*rand(No_MobileNode,2);
Distance = zeros(No_Anchors,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : No_Anchors
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the Scenario

figure('Name','Mobile Node:100, Anchor Node:4, N/W Size:300','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements

Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors,No_MobileNode)-1/2);
No_Iteration = 5;

MobileLoc_Est = NetworkSize3*rand(No_MobileNode,2);     % Initial Guess (Random Locatio)

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)-Anchor_Node_Loc(:,1))./Dist_Esti ...      % x-coordinate
            (MobileLoc_Est(m,2)-Anchor_Node_Loc(:,2))./Dist_Esti];                      % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',8,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location',...
    'Location','NorthEast')

% Compute the Root Mean Squred Error
Err3 = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/100;
axis([-0.1 1.1 -0.1 1.1]*NetworkSize3)

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, Anchor Node:4, N/W Size:300','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%% Setting Parameters

No_Anchors = 4;                 % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
NetworkSize4 = 400;             % We consider a 400by400 area that the mobile can wander

Anchor_Node_Loc   = [0   0;     % set the Anchor at 4 vertices of the region
    NetworkSize4           0;
    0           NetworkSize4;
    NetworkSize4 NetworkSize4];

Mobile_Node_Loc  = NetworkSize4*rand(No_MobileNode,2);
Distance = zeros(No_Anchors,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : No_Anchors
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the Scenario

figure('Name','Mobile Node:100, Anchor Node:4, N/W Size:400','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements

Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors,No_MobileNode)-1/2);
No_Iteration = 5;

MobileLoc_Est = NetworkSize4*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)-Anchor_Node_Loc(:,1))./Dist_Esti ...          % x-coordinate
            (MobileLoc_Est(m,2)-Anchor_Node_Loc(:,2))./Dist_Esti];                          % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',8,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location',...
    'Location','NorthEast')

% Compute the Root Mean Squred Error
Err4 = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/100;
axis([-0.1 1.1 -0.1 1.1]*NetworkSize4)


f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, Anchor Node:4, N/W Size:400','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%% Setting Parameters

No_Anchors = 4;                 % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
NetworkSize5 = 500;             % We consider a 500by500 area that the mobile can wander

Anchor_Node_Loc   = [0     0;   % set the Anchor at 4 vertices of the region
    NetworkSize5           0;
    0           NetworkSize5;
    NetworkSize5 NetworkSize5];

Mobile_Node_Loc  = NetworkSize5*rand(No_MobileNode,2);
Distance = zeros(No_Anchors,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : No_Anchors
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the Scenario

figure('Name','Mobile Node:100, Anchor Node:4, N/W Size:500','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements

Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors,No_MobileNode)-1/2);
No_Iteration = 5;

MobileLoc_Est = NetworkSize5*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)-Anchor_Node_Loc(:,1))./Dist_Esti ...      % x-coordinate
            (MobileLoc_Est(m,2)-Anchor_Node_Loc(:,2))./Dist_Esti];                      % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',8,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location',...
    'Location','NorthEast')

% Compute the Root Mean Squred Error
Err5 = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/100;
axis([-0.1 1.1 -0.1 1.1]*NetworkSize5)


f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, Anchor Node:4, N/W Size:500','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%% Number Of Nodes Vs Localization NRMSE[m] (Graph)

Err=[Err5 Err4 Err3 Err2 Err1];
Network_Size=[NetworkSize1 NetworkSize2 NetworkSize3 NetworkSize4 NetworkSize5];

RALx=[100 200 300 400 500];
RALy=[1 .82 .75 .7 .68];

DVx=[100 200 300 400 500];
DVy=[1.3 .93 .87 .85 .79];

figure('Name','Number Of Nodes Vs Localization NRMSE[m]','NumberTitle','off');
xlabel('Number Of Nodes');
ylabel('Localization NRMSE[m]');
hold on
grid on
plot(Network_Size, Err, 'R-^');
plot(RALx,RALy,'B-d');
plot(DVx,DVy,'G-s');
legend('Proposed Algorithm','RAL','DV-HOP');
hold off

%% Setting Parameters

No_Anchors2 = 2;                % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
Network_Size = 100;             % We consider a 100by100 area that the mobile can wander

Anchor_Node_Loc   = [0     0;   % Set the Anchor at 2 vertices of the region
    Network_Size       Network_Size];
Mobile_Node_Loc  = Network_Size*rand(No_MobileNode,2);
Distance = zeros(No_Anchors2,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : No_Anchors2
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the scenario
figure('Name','Mobile Node:100, N/W Size:100, Anchor Node:2','NumberTitle','off');
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements
Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors2,No_MobileNode)-1/2);
No_Iteration = 5;

% Initial Guess (random locatio)
MobileLoc_Est = Network_Size*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors2,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)- Anchor_Node_Loc(:,1))./Dist_Esti ...         % X-coordinate
            (MobileLoc_Est(m,2)- Anchor_Node_Loc(:,2))./Dist_Esti];                         % Y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',10,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location','Location','NorthEast')

% Compute the Root Mean Squred Error
Erra = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/1000;
axis([-0.1 1.1 -0.1 1.1]*Network_Size)

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, N/W Size:100, Anchor Node:2','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%% Setting Parameters

N2 = 4;                         % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes
Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
Network_Size = 100;             % We Consider a 100by100 area that the mobile can wander

Anchor_Node_Loc   = [0  0;      % Set the Anchor at 2 vertices of the region
    Network_Size/2          0;
    0           Network_Size;
    Network_Size  Network_Size/2];

Mobile_Node_Loc  = Network_Size*rand(No_MobileNode,2);
Distance = zeros(N2,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : N2
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the scenario
figure('Name','Mobile Node:100, N/W Size:100, Anchor Node:4','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements
Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(N2,No_MobileNode)-1/2);
No_Iteration = 5;

MobileLoc_Est = Network_Size*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),N2,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)- Anchor_Node_Loc(:,1))./Dist_Esti ...         % x-coordinate
            (MobileLoc_Est(m,2)- Anchor_Node_Loc(:,2))./Dist_Esti];                         % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',10,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location','Location','NorthEast')

% Compute the Root Mean Squred Error
Errb = mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2)))/100;
axis([-0.1 1.1 -0.1 1.1]*Network_Size)

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, N/W Size:100, Anchor Node:4','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%%
N3 = 6;                         % Number of Anchors
No_MobileNode = 100;            % Number of Mobile Nodes

Dist_Measure_ErrRatio = 0.1;    % Distance dependent err (standard deviation of the noise normalized to distance)
Network_Size = 100;             % We consider a 100by100 area that the mobile can wander

Anchor_Node_Loc   = [0   0;     % Set the anchor at 6 Vertices of the region
    Network_Size/2          0;
    Network_Size/2           Network_Size;
    Network_Size/2 Network_Size/2;
    0 Network_Size;
    Network_Size Network_Size];

Mobile_Node_Loc  = Network_Size*rand(No_MobileNode,2);
Distance = zeros(N3,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : N3
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the Scenario

figure('Name','Mobile Node:100, N/W Size:100, Anchor Node:6','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements
Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(N3,No_MobileNode)-1/2);
No_Iteration = 5;

MobileLoc_Est = Network_Size*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),N3,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)- Anchor_Node_Loc(:,1))./Dist_Esti ...         % x-coordinate
            (MobileLoc_Est(m,2)- Anchor_Node_Loc(:,2))./Dist_Esti];                         % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',10,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location','Location','NorthEast')

% Compute the Root Mean Squred Error
Errc = mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2)))/100;
axis([-0.1 1.1 -0.1 1.1]*Network_Size)

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, N/W Size:100, Anchor Node:6','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%%
N4 = 8;                         % Number of anchors
No_MobileNode = 100;            % Number of mobile nodes
Dist_Measure_ErrRatio = 0.1;    % distance dependent err (standard deviation of the noise normalized to distance)
Network_Size = 100;             % We consider a 100by100 area that the mobile can wander

Anchor_Node_Loc   = [0   0;     % Set the anchor at 8 vertices of the region
    Network_Size          0;
    Network_Size/2           Network_Size;
    Network_Size/2 Network_Size/2;
    0 Network_Size;
    Network_Size Network_Size;
    Network_Size/2 Network_Size;
    Network_Size/2 0];

Mobile_Node_Loc  = Network_Size*rand(No_MobileNode,2);
Distance = zeros(N4,No_MobileNode);

for m = 1 : No_MobileNode
    for n = 1 : N4
        Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
            (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
    end
end

% Plot the scenario
figure('Name','Mobile Node:100, N/W Size:100, Anchor Node:8','NumberTitle','off');clf
plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);

% Noisy Measurements
Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(N4,No_MobileNode)-1/2);
No_Iteration = 5;

MobileLoc_Est = Network_Size*rand(No_MobileNode,2);

for m = 1 : No_MobileNode
    for i = 1 : No_Iteration
        Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),N4,1)).^2 , 2));
        Distance_Drv   = [(MobileLoc_Est(m,1)- Anchor_Node_Loc(:,1))./Dist_Esti ... % x-coordinate
            (MobileLoc_Est(m,2)- Anchor_Node_Loc(:,2))./Dist_Esti];   % y-coordinate
        Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
        MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
    end
end

plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',10,'lineWidth',2);
xlabel('X axis');
ylabel('Y axis');
legend('Anchor locations','Mobile true location','Mobile estimated location','Location','NorthEast')

% Compute the Root Mean Squred Error
Errd = mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2)))/100;
axis([-0.1 1.1 -0.1 1.1]*Network_Size)
hold off

f=figure('Name','Distance Between All the Nodes and  Anchor Nodes in the network','NumberTitle','off');
t = uitable('Parent',f,'Data',Distance);
title('Mobile Node:100, N/W Size:100, Anchor Node:8','fontsize',12,'fontname','Times New Roman','color','Black');
axis off

%%

ERR=[Erra Errb Errc Errd];
No_Anchors=[No_Anchors2 N2 N3 N4];
RALa=[2 4 6 8];
RALb=[.86 .82 .8 .81];
DVa=[2 4 6 8];
DVb=[.98 .99 .93 .94];

figure('Name','Anchor Percentage Vs Localization NRMSE[m]','Numbertitle','off');
hold on
grid
plot(No_Anchors,ERR,'R-^');
plot(RALa,RALb,'B-d');
plot(DVa,DVb,'g-s');
xlabel('Anchor Percentage');
ylabel('Localization NRMSE[m]');
legend('Proposed Algorithm','RAL','DV-HOP');
hold off

%%

No_MobileNode = 100;                % Number of mobile nodes
No_Anchors = 4;                     % Number of anchors

Dist_Measure_ErrRatio = 0.1;        % Distance dependent err (standard deviation of the noise normalized to distance)

network=400:-25:50
for iii= 1:length(network)
    Network_Size = network(iii);
    Anchor_Node_Loc   = [0   0;     % set the anchor at 4 vertices of the region
        Network_Size           0;
        0           Network_Size;
        Network_Size Network_Size];
    
    Mobile_Node_Loc  = Network_Size*rand(No_MobileNode,2);
    Distance = zeros(No_Anchors,No_MobileNode);
    
    for m = 1 : No_MobileNode
        for n = 1 : No_Anchors
            Distance(n,m) = sqrt( (Anchor_Node_Loc(n,1)-Mobile_Node_Loc(m,1)).^2 + ...
                (Anchor_Node_Loc(n,2)-Mobile_Node_Loc(m,2)).^2  );
        end
    end
    
    % Plot the scenario
    figure(100);    clf
    plot(Anchor_Node_Loc(:,1),Anchor_Node_Loc(:,2),'bo','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
    grid on
    hold on
    plot(Mobile_Node_Loc(:,1),Mobile_Node_Loc(:,2),'ks','MarkerSize',8,'lineWidth',2);
    
    % noisy measurements
    Distance_Noisy = Distance + Distance.*Dist_Measure_ErrRatio.*(rand(No_Anchors,No_MobileNode)-1/2);
    No_Iteration = 5;
    MobileLoc_Est = Network_Size*rand(No_MobileNode,2);
    
    for m = 1 : No_MobileNode
        for j = 1 : No_Iteration
            Dist_Esti   = sqrt(sum( (Anchor_Node_Loc - repmat(MobileLoc_Est(m,:),No_Anchors,1)).^2 , 2));
            Distance_Drv   = [(MobileLoc_Est(m,1)-Anchor_Node_Loc(:,1))./Dist_Esti ... % x-coordinate
                (MobileLoc_Est(m,2)-Anchor_Node_Loc(:,2))./Dist_Esti];   % y-coordinate
            Delta = - (Distance_Drv.'*Distance_Drv)^-1*Distance_Drv.' * (Dist_Esti - Distance_Noisy(:,m));
            MobileLoc_Est(m,:) = MobileLoc_Est(m,:) + Delta.';
        end
    end
    
    plot(MobileLoc_Est(:,1),MobileLoc_Est(:,2),'m^','MarkerSize',8,'lineWidth',2);
    xlabel('X axis');
    ylabel('Y axis');
    legend('Anchor locations','Mobile true location','Mobile estimated location',...
        'Location','NorthEast')
    
    % Compute the Root Mean Squred Error
    Err(iii) = (mean(sqrt(sum((MobileLoc_Est-Mobile_Node_Loc).^2))))/100;
    axis([-0.1 1.1 -0.1 1.1]*Network_Size)
end

%% Localization NRMSE[m] Vs CDF

R=75:-2.3:1
for j=1:length(R)
    x=5
    r=5
    No_Anchors=4
    S=3000
    lambda=No_Anchors/S
    theta=acos(x/mean(2*R))
    theta1=acos(x/(2*r))
    a=(theta-(sin(2*theta)/2))
    b=(theta1-(sin(2*theta1)/2))
    f=(R.^2)*(a+b)
    fz=exp(-(lambda*f))
    FZ=fliplr(fz)
end

for k=1:15
    ERR=Err(1:k);
    fz1=FZ(1:k)
end

figure('Name','Localization NRMSE[m] Vs CDF','NumberTitle','off');
grid on
hold on
xlabel('Localization NRMSE[m]');
ylabel('CDF');

DVcdfx=[0 .04 .08 .15 .2 .23 .3 .35 .4 .43 .5 .55 .6 .64 .7 .75 .8 ];
DVcdfy=[0  0 .01 .01 .03 .05 .075 .1 .15 .2 .3 .4 .5  .59  .9 1 1 ];

RALcdfx=[0 .05 .1 .15 .2 .25 .3 .35 .4 .44 .5 .55 .6 .65 .7 .75 .8 .85 .9];
RALcdfy=[0 0 .02 .04 .05 .1 .15 .25 .35 .45 .6 .75 .82  .9 .95 1 1 1 1];

plot(ERR,fz1,'r^--','MarkerSize',8,'lineWidth',2);
plot(RALcdfx,RALcdfy,'bd-','MarkerSize',8,'lineWidth',2);
plot(DVcdfx,DVcdfy,'gs-','MarkerSize',8,'lineWidth',2);
legend('Proposed Algorithm', ' RAL','DV-HOP','Location','SouthEast');
hold off

%%