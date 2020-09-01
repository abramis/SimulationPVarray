clc
clear all;
    %% Problem definition
 
    CostFunction =@(Irow1,Irow2,Irow3,maxRowC,panelMax) SphereArxiko(Irow1,Irow2,Irow3,maxRowC); %Cost Function
    nVar = 1;     % Number of Unknown (Decision) Variables
    VarSize = [1 nVar];      % Matrix Size of Decision Variables

    VarMin = 0.1; % Lower Bound of Decision Variable
    VarMax = 9; % Upeer Bound of Decision Variables

    %% Parameters of PSO

    MaxIt = 100;  % Maximum number of iterations

    nPop = 9;  % Population size (Swarm Size)

    w = 1;  % Inertsia Coefficient
    wdamp = 0.99; %Dampin Ration of Inertia Weight
    c1 = 1.2;  % Personal Acceleration Coefficient
    c2 = 1.8;  % Social Acceleration Coefficient
    
    %The Flag for Showing Iteration information
    ShowIterInfo = true; 
    
    %Define the Bounds of Velocities
   MaxVelocity = (VarMax-VarMin)*0.5;
   MinVelocity = -MaxVelocity;
    
    %% Initialization
    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Current=[];
    empty_particle.Coordinates=[];
    empty_particle.Irradiation=[];
    empty_particle.Voltage=[];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];
    empty_particle.Best.Current=[];
 

    %Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    %Initialize Global Best
    GlobalBest.Cost = inf;

%% Initialize Irradiations in the Array
Irr(1,1)=600;
Irr(1,2)=800;
Irr(1,3)=900;

Irr(2,1)=400;
Irr(2,2)=800;
Irr(2,3)=900;

Irr(3,1)=200;
Irr(3,2)=800;
Irr(3,3)=900;

% Run the Simulink
sim('tepeRarx.slx');


%% Plot P-I 
P=Vtot.*Itot;
yplot2=P;
yplot2(yplot2<0)=nan;
figure;
plot(Vtot,yplot2);
xlabel('Voltage(V)');
ylabel('Power(W)');

%%Find Vmp, Imp Total
PmaxFromTot=max(P);
[PowerMax,indx]=max(P);
Vmp=Vtot(indx);
Imp=Itot(indx);
Pmax=Vmp*Imp;
%%Plot I-V
%%Plot I-V
yplot3=Itot;
yplot3(yplot3<0)=nan;
figure;
plot(Vtot,yplot3);
xlabel('Voltage');
ylabel('Current');

% Find maxV, maxI for every row and total
maxVtot=max(Vtot);
maxItot=max(Itot);
% maxI1=max(I1);
% maxI2=max(I2);
% maxI3=max(I3);

% Find Imp for every row
% I1mp=I1(indx);
% I2mp=I2(indx);
% I3mp=I3(indx);

% Compute Imp, Vmp for every module

%%11
maxPmod11Res=mod11(:,1).*mod11(:,2);
[maxPmod11,indMod11]=max(maxPmod11Res);
Imp(1,1) = mod11(indMod11,2);
Vmp(1,1) = mod11(indMod11,1);
P(1,1) = Imp(1,1)*Vmp(1,1);

%%12
maxPmod12Res=mod12(:,1).*mod12(:,2);
[maxPmod12,indMod12]=max(maxPmod12Res);
Imp(1,2) = mod12(indMod12,2);
Vmp(1,2) = mod12(indMod12,1);
P(1,2) = Imp(1,2)*Vmp(1,2);

%13
maxPmod13Res=mod13(:,1).*mod13(:,2);
[maxPmod13,indMod13]=max(maxPmod13Res);
Imp(1,3) = mod13(indMod13,2);
Vmp(1,3) = mod13(indMod13,1);
P(1,3) = Imp(1,3)*Vmp(1,3);

%21
maxPmod21Res=mod21(:,1).*mod21(:,2);
[maxPmod21,indMod21]=max(maxPmod21Res);
Imp(2,1) = mod21(indMod21,2);
Vmp(2,1) = mod21(indMod21,1);
P(2,1) = Imp(2,1)*Vmp(2,1);

%%22
maxPmod22Res=mod22(:,1).*mod22(:,2);
[maxPmod22,indMod22]=max(maxPmod22Res);
Imp(2,2) = mod22(indMod22,2);
Vmp(2,2) = mod22(indMod22,1);
P(2,2) = Imp(2,2)*Vmp(2,2);

%%23
maxPmod23Res=mod23(:,1).*mod23(:,2);
[maxPmod23,indMod23]=max(maxPmod23Res);
Imp(2,3) = mod23(indMod23,2);
Vmp(2,3) = mod23(indMod23,1);
P(2,3) = Imp(2,3)*Vmp(2,3);

%%31
maxPmod31Res=mod31(:,1).*mod31(:,2);
[maxPmod31,indMod31]=max(maxPmod31Res);
Imp(3,1) = mod31(indMod31,2);
Vmp(3,1) = mod31(indMod31,1);
P(3,1) = Imp(3,1)*Vmp(3,1);


%%32
maxPmod32Res=mod32(:,1).*mod32(:,2);
[maxPmod32,indMod32]=max(maxPmod32Res);
Imp(3,2) = mod32(indMod32,2);
Vmp(3,2) = mod32(indMod32,1);
P(3,2) = Imp(3,2)*Vmp(3,2);

%%33
maxPmod33Res=mod33(:,1).*mod33(:,2);
[maxPmod33,indMod33]=max(maxPmod33Res);
Imp(3,3) = mod33(indMod33,2);
Vmp(3,3) = mod33(indMod33,1);
P(3,3) = Imp(3,3)*Vmp(3,3);

Icolumn(1)= Imp(1,1)+Imp(2,1)+Imp(3,1);
Icolumn(2)= Imp(1,2)+Imp(2,2)+Imp(3,2);
Icolumn(3)= Imp(1,3)+Imp(2,3)+Imp(3,3);

Ibegin(1)=Icolumn(1);
Ibegin(2)=Icolumn(2);
Ibegin(3)=Icolumn(3);
%Find max row current
maxRowCurrent=Icolumn(1);
for i=2:3
    if(Icolumn(i)>maxRowCurrent)
        maxRowCurrent=Icolumn(i);
    end
end
PmaxNoBypassing=Vmp(1,1)*Icolumn(1)+Vmp(1,2)*Icolumn(2)+Vmp(1,3)*Icolumn(3);
    
    %% Initialize Population Members

        %Currents
        particle(1).Position = 1;  
        particle(2).Position = 2;  
        particle(3).Position = 3;   
        particle(4).Position = 4;
        particle(5).Position = 5;
        particle(6).Position = 6;
        particle(7).Position = 7;
        particle(8).Position = 8;
        particle(9).Position = 9;
        
        particle(1).Coordinates = 11;
        particle(2).Coordinates = 21;
        particle(3).Coordinates = 31;
        particle(4).Coordinates = 12;
        particle(5).Coordinates = 22;
        particle(6).Coordinates = 32;
        particle(7).Coordinates = 13;
        particle(8).Coordinates = 23;
        particle(9).Coordinates = 33;

        %Irradiation
        particle(1).Irradiation = Irr(1,1);
        particle(2).Irradiation = Irr(2,1);
        particle(3).Irradiation = Irr(3,1);
        particle(4).Irradiation = Irr(1,2);
        particle(5).Irradiation = Irr(2,2);
        particle(6).Irradiation = Irr(3,2);
        particle(7).Irradiation = Irr(1,3);
        particle(8).Irradiation = Irr(2,3);
        particle(9).Irradiation = Irr(3,3);
        
        %Currents
        particle(1).Current = Imp(1,1);
        particle(2).Current = Imp(2,1);
        particle(3).Current = Imp(3,1);
        particle(4).Current = Imp(1,2);
        particle(5).Current = Imp(2,2);
        particle(6).Current = Imp(3,2);
        particle(7).Current = Imp(1,3);
        particle(8).Current = Imp(2,3);
        particle(9).Current = Imp(3,3);
        
        
        %Voltages
        particle(1).Voltage = Vmp(1,1);
        particle(2).Voltage = Vmp(2,1);
        particle(3).Voltage = Vmp(3,1);
        particle(4).Voltage = Vmp(1,2);
        particle(5).Voltage = Vmp(2,2);
        particle(6).Voltage = Vmp(3,2);
        particle(7).Voltage = Vmp(1,3);
        particle(8).Voltage = Vmp(2,3);
        particle(9).Voltage = Vmp(3,3);
        

        
        for i =1:nPop
            %Initialize Velocity
            particle(i).Velocity = 1+round(rand()*8);

            % Evaluation
            particle(i).Cost = CostFunction(Icolumn(1),Icolumn(2),Icolumn(3),maxRowCurrent);
        
   
        
            %Update the Personal Best
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost= particle(i).Cost;
            particle(i).Best.Current=particle(i).Current;
            particle(i).Best.Coordinates=particle(i).Coordinates;
            particle(i).Best.Irradiation=particle(i).Irradiation;

            %Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end


        end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt,1);

    %% Main Loop of PSO

    for it =1 : MaxIt

        for i=1:nPop

            %Update Velocity
            particle(i).Velocity = w*particle(i).Velocity...
                + c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position)...
                + c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);

            %Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);            
            
            %Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            %Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);
            

            
            % Evaluation
            particle(i).Cost = CostFunction(Icolumn(1),Icolumn(2),Icolumn(3),maxRowCurrent);
        

            %Update Personal Best
            if(particle(i).Cost < particle(i).Best.Cost)
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;
                particle(i).Best.Current=particle(i).Current;
                particle(i).Best.Coordinates=particle(i).Coordinates;
                particle(i).Best.Irradiation=particle(i).Irradiation;
                

                %Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end

            end



        end

        %Store the Best Cost Value
        BestCosts(it) = GlobalBest.Cost;
        
        %Display Iteration Information
        if(ShowIterInfo)
            disp(['Iteration' num2str(it) ': Best Cost =' num2str(BestCosts(it))]);
        end
        %Damping Inertia Coefficient
        w = w*wdamp;
        
        for inner = 1 : nPop
            k(inner)=particle(inner).Position;
        end
        [value,index]=sort(k);
        for i=1:nPop
            particle(i).Position = value(i);        
            tempIrradiation(i)=particle(index(i)).Irradiation;
            tempCoordinates(i)=particle(index(i)).Coordinates;
            tempCurrent(i)=particle(index(i)).Current;
        end
        
        for i=1:nPop
            particle(i).Current = tempCurrent(i);
            particle(i).Coordinates = tempCoordinates(i);
            particle(i).Irradiation = tempIrradiation(i);
        end
        
        
        Icolumn(1)= particle(1).Current+particle(2).Current+particle(3).Current;
        Icolumn(2)= particle(4).Current+particle(5).Current+particle(6).Current;
        Icolumn(3)= particle(7).Current+particle(8).Current+particle(9).Current;

            maxRowCurrent=Icolumn(1);
        for i=2:3
            if(Icolumn(i)>maxRowCurrent)
                maxRowCurrent=Icolumn(i);
            end
        end
    end
    
    
    IcolumnBest(1)=particle(1).Best.Current+particle(2).Best.Current+particle(3).Best.Current;
    IcolumnBest(2)=particle(4).Best.Current+particle(5).Best.Current+particle(6).Best.Current;
    IcolumnBest(3)=particle(7).Best.Current+particle(8).Best.Current+particle(9).Best.Current;
    
    
   
        Irr(1,1)=particle(1).Best.Irradiation;
        Irr(2,1)=particle(2).Best.Irradiation;
        Irr(3,1)=particle(3).Best.Irradiation;
        Irr(1,2)=particle(4).Best.Irradiation;
        Irr(2,2)=particle(5).Best.Irradiation;
        Irr(3,2)=particle(6).Best.Irradiation;
        Irr(1,3)=particle(7).Best.Irradiation;
        Irr(2,3)=particle(8).Best.Irradiation;
        Irr(3,3)=particle(9).Best.Irradiation;
        
    
%Run the Simulink
sim('tepeRarx.slx');
    
    yplot1=Itot;
    yplot1(yplot1<0)=nan;
    figure(3);
    plot(Vtot,yplot1,'m:','LineWidth',2);
    xlabel('Voltage(V)');
    ylabel('Current(A)');
    hold on
    Pfinal=Vtot.*Itot;
    yplot=Pfinal;
    yplot(yplot<0)=nan;
    [Peak, PeakIdx] = max(yplot);
    figure(4);
    plot(Vtot,yplot,'m:',Vtot(PeakIdx),yplot(PeakIdx),'x','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',8);
    xlabel('Voltage(V)');
    ylabel('Power(W)');
    text(Vtot(PeakIdx), Peak, sprintf('PSO-GP=%6.1f W', Peak))
    hold on
    out.pop=particle;
    out.BestSol=GlobalBest;
    out.BestCosts=BestCosts;
    