clc
clear all;
    %% Problem definition
 
    CostFunction =@Sphere; %Cost Function

    %% Parameters of PSO

    MaxIt = 500;  % Maximum number of iterations

    %The Flag for Showing Iteration information
    ShowIterInfo = true; 
    numberOfModules=81;
    numberOfGenes=100;
    
    %% Initialization
    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Current=[];
    empty_particle.Coordinates=[];
    empty_particle.Irradiation=[];
    empty_particle.Voltage=[];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];
    empty_particle.Best.Current=[];
    
    empty_gen.Cost=[];
    empty_gen.Number=[];

    %Create Population Array
    particle = repmat(empty_particle, numberOfModules, 1);

    %Initialize Global Best
    GlobalBest.Cost = inf;

%% Initialize Irradiations in the Array
Irr(1)=600;Irr(2)=900;Irr(3)=600;
Irr(10)=600;Irr(11)=900;Irr(12)=600;
Irr(19)=600;Irr(20)=900;Irr(21)=600;

Irr(28)=400;Irr(29)=400;Irr(30)=400;
Irr(37)=200;Irr(38)=200;Irr(39)=200;
Irr(46)=400;Irr(47)=400;Irr(48)=400;

Irr(55)=200;Irr(56)=200;Irr(57)=900;
Irr(64)=200;Irr(65)=800;Irr(66)=200;
Irr(73)=800;Irr(74)=200;Irr(75)=200;

Irr(4)=800;Irr(13)=800;Irr(22)=200;Irr(31)=800;Irr(40)=800;Irr(49)=800;Irr(58)=800;Irr(67)=800;Irr(76)=800;

Irr(5)=900;Irr(14)=900;Irr(23)=900;Irr(32)=900;Irr(41)=900;Irr(50)=800;Irr(59)=900;Irr(68)=900;Irr(77)=900;
Irr(6)=900;Irr(15)=900;Irr(24)=600;Irr(33)=900;Irr(42)=900;Irr(51)=900;Irr(60)=900;Irr(69)=900;Irr(78)=600;
Irr(7)=400;Irr(16)=400;Irr(25)=900;Irr(34)=900;Irr(43)=900;Irr(52)=600;Irr(61)=900;Irr(70)=900;Irr(79)=900;
Irr(8)=400;Irr(17)=400;Irr(26)=900;Irr(35)=600;Irr(44)=900;Irr(53)=900;Irr(62)=900;Irr(71)=900;Irr(80)=900;
Irr(9)=900;Irr(18)=900;Irr(27)=900;Irr(36)=600;Irr(45)=900;Irr(54)=900;Irr(63)=400;Irr(72)=900;Irr(81)=900;

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

% Compute Imp, Vmp for every module

for i=1:numberOfModules
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
        for i=1:numberOfModules
            particle(i).Position = i;  
        end
        temp2=11;
        for i=1:numberOfModules
            particle(i).Coordinates = temp2;
            temp2=temp2+10;
            if mod(i,9)==0
                temp2=10+rdivide(i,9)+1;
            end
        end

        %Irradiation %Currents %Voltages
        temp3=1;
        j=1;
        for i=1:numberOfModules
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
        
        gen(1).Number=[particle(:).Position]; 
        % Evaluation
        gen(1).Cost = CostFunction(Icolumn,maxRowCurrent);
        
        BestCost=gen(1).Cost;
        BestGene=gen(1).Number;
        
        %% Create Random Genes
 
    for it=1:numberOfGenes
        % Initialize Position
        for j=1:numberOfModules
            particle(j).Position=j;
        end
        
        for j=1:numberOfModules
            particle(j).Coordinates = temp2;
            temp2=temp2+10;
            if mod(j,9)==0
                temp2=10+rdivide(j,9)+1;
            end
        end

        %Irradiation %Currents %Voltages
        temp3=1;
        j=1;
        for i=1:numberOfModules
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
        
        %% Make the computation
            idx=randperm(numberOfModules);
            final(it,:)=idx;
            
            for j=1:numberOfModules
                tempIrradiation(j)=particle(idx(j)).Irradiation;
                tempCurrent(j)=particle(idx(j)).Current;
                tempCoordinates(j)=particle(idx(j)).Coordinates;
                particle(j).Position=idx(j);
            end
            
            gen(it).Number=[particle(:).Position];
            for j=1:numberOfModules
                particle(j).Irradiation=tempIrradiation(j);
                particle(j).Current=tempCurrent(j);
                particle(j).Coordinates=tempCoordinates(j);
            end
            
          %Evaluation
          %Find row columns
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
         
         %Compute Cost
         gen(it).Cost = CostFunction(Icolumn,maxRowCurrent);
         
         if(gen(it).Cost < BestCost)
            BestCost=gen(it).Cost;
            BestGene=gen(it).Number;
            j=1;
            for i=1:9
                Icolumn(i)= particle(j).Current+particle(j+1).Current+particle(j+2).Current...
                    +particle(j+3).Current+particle(j+4).Current+particle(j+5).Current...
                    +particle(j+6).Current+particle(j+7).Current+particle(j+8).Current;
                j=j+9;
            end
            for h=1:numberOfModules
                particle(h).Best.Irradiation=particle(h).Irradiation;
            end
             
         end
         
        
    end
     
     
     %% Mutation, Crossover Rand
  for iter=1:MaxIt     
     %Sort the genes
         [value,index]=sort([gen(:).Cost]);
         
         for i=1:numberOfGenes
             tempGen(i,:)=gen(index(i)).Number;
         end
         
         %Hold the 40 first Genes
         %Create the Mutations
         rowNumber=41;
         for i=1:40          
                 randNum1=ceil(81*rand()); %First Index Number

                 randNum2=ceil(81*rand()); %Second Index Number
                 while(randNum2==randNum1)
                    randNum2=ceil(81*rand());
                 end
                 
                 %Swap values
                 temporalGen=tempGen(i,:);
                 tempNumber=temporalGen(randNum1);
                 temporalGen(randNum1)=temporalGen(randNum2);
                 temporalGen(randNum2)=tempNumber;
                 
                 tempGen(rowNumber,:)=temporalGen;
                 rowNumber=rowNumber+1;
                 
         end
         
         %% Create random genes
         for i=1:20
            idxNew=randperm(numberOfModules);
            finalNew(i,:)=idxNew;
            
            tempGen(rowNumber,:)=finalNew(i,:);
 
            rowNumber=rowNumber+1;
         end
         
         %% Evaluation
         
    for it=1:numberOfGenes
         %Initialize Position
        for v=1:numberOfModules
                particle(v).Position=v;
        end
        

        %Initialize Coordinates
        for i=1:numberOfModules
            particle(i).Coordinates = temp2;
            temp2=temp2+10;
            if mod(i,9)==0
                temp2=10+rdivide(i,9)+1;
            end
        end

        %Irradiation %Currents %Voltages
        temp3=1;
        j=1;
        for i=1:numberOfModules
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
         
         
         
         
            for j=1:numberOfModules
                tempPosition(j)=tempGen(it,j);
            end
            
            for j=1:numberOfModules
                tempIrradiation(j)=particle(tempPosition(j)).Irradiation;
                tempCurrent(j)=particle(tempPosition(j)).Current;
                tempCoordinates(j)=particle(tempPosition(j)).Coordinates;
                particle(j).Position=tempPosition(j);
            end
            
            gen(it).Number=[particle(:).Position];
            for j=1:numberOfModules
                particle(j).Irradiation=tempIrradiation(j);
                particle(j).Current=tempCurrent(j);
                particle(j).Coordinates=tempCoordinates(j);
            end
            
          %Evaluation
          %Find row columns
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
         
         %Compute Cost
         gen(it).Cost = CostFunction(Icolumn,maxRowCurrent);
         
         if(gen(it).Cost < BestCost)
            BestCost=gen(it).Cost;
            BestGene=gen(it).Number;
            j=1;
            for i=1:9
                IcolumnBest(i)= particle(j).Current+particle(j+1).Current+particle(j+2).Current...
                    +particle(j+3).Current+particle(j+4).Current+particle(j+5).Current...
                    +particle(j+6).Current+particle(j+7).Current+particle(j+8).Current;
                j=j+9;
            end
            for h=1:numberOfModules
                particle(h).Best.Irradiation=particle(h).Irradiation;
            end
             
         end
    end
 

          [value1,index1]=sort([gen(:).Cost]);
  end    
         %% ShowResults
       j=1; temp7=1;
    for i=1:numberOfModules
        Irr(j)=particle(i).Best.Irradiation;
        aktin(i)=particle(i).Best.Irradiation;
        j=j+9;
        if mod(i,9)==0
            temp7=temp7+1;
            j=temp7;
        end
    end     
    
aktinPinakas=reshape(aktin.',9,9).';
aktinPinakasFinal=flipud(aktinPinakas);
IcolumnFinal=flipud(reshape(IcolumnBest.',1,9).');  
%Run the Simulink
sim('power_PVArray_PartialShading9x9new.slx');
    

    yplot1=Itot;
    yplot1(yplot1<0)=nan;
    figure(3);
    plot(Vtot,yplot1,'b--','LineWidth',2);
    xlabel('Voltage(V)');
    ylabel('Current(A)');
    hold on
    Pfinal=Vtot.*Itot;
    yplot=Pfinal;
    yplot(yplot<0)=nan;
    [Peak, PeakIdx] = max(yplot);
    figure(4);
    plot(Vtot,yplot,'b--',Vtot(PeakIdx),yplot(PeakIdx),'+','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',8);
    xlabel('Voltage(V)');
    ylabel('Power(W)');
    text(Vtot(PeakIdx)+10, Peak-1500, sprintf('GA-GP=%6.1f W', Peak))
    hold on
      