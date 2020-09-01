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
Irr(1)=200;Irr(2)=200;Irr(3)=200;
Irr(10)=200;Irr(11)=200;Irr(12)=200;
Irr(19)=200;Irr(20)=200;Irr(21)=200;

Irr(28)=200;Irr(29)=200;Irr(30)=200;
Irr(37)=200;Irr(38)=200;Irr(39)=200;
Irr(46)=200;Irr(47)=200;Irr(48)=200;

Irr(55)=200;Irr(56)=200;Irr(57)=200;
Irr(64)=200;Irr(65)=200;Irr(66)=200;
Irr(73)=200;Irr(74)=200;Irr(75)=200;

Irr(4)=600;Irr(13)=600;Irr(22)=600;Irr(31)=600;Irr(40)=600;Irr(49)=600;Irr(58)=600;Irr(67)=600;Irr(76)=600;

Irr(5)=600;Irr(14)=600;Irr(23)=600;Irr(32)=600;Irr(41)=600;Irr(50)=600;Irr(59)=600;Irr(68)=600;Irr(77)=600;
Irr(6)=600;Irr(15)=600;Irr(24)=600;Irr(33)=600;Irr(42)=600;Irr(51)=600;Irr(60)=600;Irr(69)=600;Irr(78)=600;
Irr(7)=900;Irr(16)=900;Irr(25)=900;Irr(34)=900;Irr(43)=900;Irr(52)=900;Irr(61)=900;Irr(70)=900;Irr(79)=900;
Irr(8)=900;Irr(17)=900;Irr(26)=900;Irr(35)=900;Irr(44)=900;Irr(53)=900;Irr(62)=900;Irr(71)=900;Irr(80)=900;
Irr(9)=900;Irr(18)=900;Irr(27)=900;Irr(36)=900;Irr(45)=900;Irr(54)=900;Irr(63)=900;Irr(72)=900;Irr(81)=900;

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
        
         for it=1:numberOfGenes
            xGen=final(it,:);
            tempF=final;
            tempF(it,:)=[];
            tempABC=datasample(tempF,3);
            aGen=tempABC(1,:);
            bGen=tempABC(2,:);
            cGen=tempABC(3,:);
            
            for j=1:numberOfModules
                r=rand();
                F=rand()*0.5;
                flagCheck=false;
                %Check for crossover propability
                if r>0.5
                    yGen(j)=round(aGen(j)+abs(F*(bGen(j)-cGen(j))));
                    %Check for same values
                    for k=1:j
                        if(yGen(j)==yGen(k) && j~=k)
                            flagCheck=true;
                            break;
                        end
                    end
                    %Check for out of range and act if something is wrong
                    flagIncrement=true;
                    flagInternal=false;
                    while(yGen(j)>=numberOfModules+1 || flagCheck==true)
                         flagInternal=false;
                        if(yGen(j)>=numberOfModules)
                            flagIncrement=false;
                        end
                        if(yGen(j)<numberOfModules && flagIncrement==true)
                            yGen(j)=yGen(j)+1;
                        end
                        
                        if(yGen(j)>1 && flagIncrement==false)
                            yGen(j)=yGen(j)-1;
                        end
                        
                        for k=1:j
                                if(yGen(j)==yGen(k) && j~=k)
                                    flagInternal=true;
                                    flagCheck=true;
                                end
                        end
                        
                        if(flagInternal==false)
                            flagCheck=false;
                        end     
                    end  
                else
                     %IF crossover propability is below threshold
                     yGen(j)=xGen(j);
                     for k=1:j
                        if(yGen(j)==yGen(k) && j~=k)
                            flagCheck=true;
                            break;
                        end
                     end
                     flagIncrement=true;
                     flagInternal=false;
                     while(yGen(j)>=numberOfModules+1 || flagCheck==true)
                         flagInternal=false;
                         if(yGen(j)>=numberOfModules)
                            flagIncrement=false;
                         end
                         if(yGen(j)<numberOfModules && flagIncrement==true)
                            yGen(j)=yGen(j)+1;
                         end
                        
                         if(yGen(j)>1 && flagIncrement==false)
                            yGen(j)=yGen(j)-1;
                         end
                        
                         for k=1:j
                                if(yGen(j)==yGen(k) && j~=k)
                                    flagInternal=true;
                                    flagCheck=true;
                                end
                         end
                        
                         if(flagInternal==false)
                            flagCheck=false;
                         end     
                     end  
                end   
            end
            
          
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
                tempPosition(j)=yGen(j);
            end
            
            for j=1:numberOfModules
                tempIrradiation(j)=particle(tempPosition(j)).Irradiation;
                tempCurrent(j)=particle(tempPosition(j)).Current;
                tempCoordinates(j)=particle(tempPosition(j)).Coordinates;
                particle(j).Position=tempPosition(j);
            end
            
            tempNumber=[particle(:).Position];
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
         tempCost = CostFunction(Icolumn,maxRowCurrent);
         
         if(tempCost < gen(it).Cost)
            final(it,:)=yGen;
            gen(it).Cost=tempCost;
            gen(it).Number=tempNumber;
         end
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
         
         
    
         %TELOS numberOfGenes           
         end 
            
   
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
    plot(Vtot,yplot1,'r-.','LineWidth',2);
    xlabel('Voltage(V)');
    ylabel('Current(A)');
    title('I-V after Optimizations')
    hold off
    leg=legend('PSO','GA','DE');
    Pfinal=Vtot.*Itot;
    yplot=Pfinal;
    yplot(yplot<0)=nan;
    [Peak, PeakIdx] = max(yplot);
    figure(4);
    plot(Vtot,yplot,'r-.',Vtot(PeakIdx),yplot(PeakIdx),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',8);
    xlabel('Voltage(V)');
    ylabel('Power(W)');
    title('P-V after Optimizations')
    text(Vtot(PeakIdx), Peak-2000, sprintf('DE-GP=%6.1f W', Peak))
    hold off
    leg=legend('PSO','PSO-GP','GA','GA-GP','DE','DE-GP');
    leg.Location=('northwest');
    close_system('power_PVArray_PartialShading9x9new.slx');
      