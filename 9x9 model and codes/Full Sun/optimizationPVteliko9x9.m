clc
clear all;
    %% Problem definition
 
    CostFunction =@(Irow,maxRowC) Sphere(Irow,maxRowC); %Cost Function
    nVar = 1;     % Number of Unknown (Decision) Variables
    VarSize = [1 nVar];      % Matrix Size of Decision Variables

    VarMin = 0.1; % Lower Bound of Decision Variable
    VarMax = 81; % Upeer Bound of Decision Variables

    %% Parameters of PSO

    MaxIt = 500;  % Maximum number of iterations

    nPop = 81;  % Population size (Swarm Size)

    w = 1;  % Inertsia Coefficient
    wdamp = 0.99; %Dampin Ration of Inertia Weight
    c1 = 1.2;  % Personal Acceleration Coefficient
    c2 = 1.8;  % Social Acceleration Coefficient
    
    %The Flag for Showing Iteration information
    ShowIterInfo = true; 
    
    %Define the Bounds of Velocities
   MaxVelocity = (VarMax-VarMin)*0.3;
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
for i=1:nPop
    Irr(i)=900;
end

% Run the Simulink
sim('power_PVArray_PartialShading9x9new.slx');


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
yplot3=Itot;
yplot3(yplot3<0)=nan;
figure;
plot(Vtot,yplot3);
xlabel('Voltage(V)');
ylabel('Current(A)');

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

for i=1:nPop
    %sim(power_PVArray_PartialShading3x3.slx);
    temp=eval(['V' num2str(i)]);
    Pmod=temp(:,1).*temp(:,2);
    [maxPmod,indMod]=max(Pmod);
    Imp(i) = temp(indMod,2);
    Vmp(i) = temp(indMod,1);
    P(i) = Imp(i)*Vmp(i);
end
for i=1:9
    Icolumn(i)= Imp(i)+Imp(i+9)+Imp(i+18)+Imp(i+27)+Imp(i+36)+Imp(i+45)+Imp(i+54)+Imp(i+63)+Imp(i+72);
end
for i=1:9
    Ibegin(i)=Icolumn(i);
end
%Find max row current
maxRowCurrent=Icolumn(1);
for i=2:9
    if(Icolumn(i)>maxRowCurrent)
        maxRowCurrent=Icolumn(i);
    end
end
tempP=0;
for i=1:9
    tempP=tempP+Vmp(i)*Icolumn(i);
end
PmaxNoBypassing=tempP;
    
    %% Initialize Population Members

        %Currents
        for i=1:nPop
            particle(i).Position = i;  
        end
        temp2=11;
        for i=1:nPop
            particle(i).Coordinates = temp2;
            temp2=temp2+10;
            if mod(i,9)==0
                temp2=10+rdivide(i,9)+1;
            end
        end

        %Irradiation %Currents %Voltages
        temp3=1;
        j=1;
        for i=1:nPop
            temp4=9*(j-1)+temp3;
            particle(i).Irradiation = Irr(temp4);
            particle(i).Current = Imp(temp4);
            particle(i).Voltage = Vmp(temp4);
            j=j+1;
            if mod(i,9)==0
                temp3=temp3+1;
                j=1;
            end
        end
        
        for i =1:nPop
            %Initialize Velocity
            particle(i).Velocity = 1+round(rand()*8);

            % Evaluation
            particle(i).Cost = CostFunction(Icolumn,maxRowCurrent);
        
   
        
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
            particle(i).Cost = CostFunction(Icolumn,maxRowCurrent);
        

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
        j=1;
        for i=1:9
            Icolumn(i)= particle(j).Current+particle(j+1).Current+particle(j+2).Current...
                +particle(j+3).Current+particle(j+4).Current+particle(j+5).Current...
                +particle(j+6).Current+particle(j+7).Current+particle(j+8).Current;
            j=j+9;
        end

            maxRowCurrent=Icolumn(1);
        for i=2:9
            if(Icolumn(i)>maxRowCurrent)
                maxRowCurrent=Icolumn(i);
            end
        end
    end
    
    j=1;
    for i=1:9
        IcolumnBest(i)= particle(j).Best.Current+particle(j+1).Best.Current+particle(j+2).Best.Current...
            +particle(j+3).Best.Current+particle(j+4).Best.Current+particle(j+5).Best.Current...
            +particle(j+6).Best.Current+particle(j+7).Best.Current+particle(j+8).Best.Current;
        j=j+9;
    end
   
    
    for i=1:nPop
        Irr(i)=particle(i).Best.Irradiation;
    end
        
    
%Run the Simulink
sim('power_PVArray_PartialShading9x9new.slx');

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
    text(Vtot(PeakIdx)-70, Peak-700, sprintf('PSO-GP=%6.1f W', Peak))
    hold on
    out.pop=particle;
    out.BestSol=GlobalBest;
    out.BestCosts=BestCosts;
    