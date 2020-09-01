clc
clear all;
    %% Problem definition
 
    CostFunction =@SphereArxiko; %Cost Function

    %% Parameters of PSO

    MaxIt = 100;  % Maximum number of iterations

    nPop = 9;  % Population size (Swarm Size)

    %The Flag for Showing Iteration information
    ShowIterInfo = true; 
    numberOfModules=9;
    numberOfGenes=50;
    
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
    particle = repmat(empty_particle, nPop, 1);

    %Initialize Global Best
    GlobalBest.Cost = inf;

%% Initialize Irradiations in the Array
Irr(1,1)=900;
Irr(1,2)=400;
Irr(1,3)=900;

Irr(2,1)=400;
Irr(2,2)=800;
Irr(2,3)=600;

Irr(3,1)=200;
Irr(3,2)=200;
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
yplot3=Itot;
yplot3(yplot3<0)=nan;
figure;
plot(Vtot,yplot3);
xlabel('Voltage');
ylabel('Current');

% Find maxV, maxI for every row and total
maxVtot=max(Vtot);
maxItot=max(Itot);


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

        %Initialize Position
        for i=1:numberOfModules
            particle(i).Position=i;  
        end
         

        %Initialize Coordinates
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
        rows=1;
        cols=1;
        for i=1:numberOfModules
                particle(i).Irradiation= Irr(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end

        %Current
         rows=1;
        cols=1;
        for i=1:numberOfModules
                particle(i).Current= Imp(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end
        
  
        %Voltages
        rows=1;
        cols=1;
        for i=1:numberOfModules
                particle(i).Voltage= Vmp(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end
        
        gen(1).Number=[particle(:).Position]; 
        % Evaluation
        gen(1).Cost = CostFunction(Icolumn(1),Icolumn(2),Icolumn(3),maxRowCurrent);
        
        
        BestCost=gen(1).Cost;
        BestGene=gen(1).Number;
        
        %% Create Random Genes
 
     for i=1:numberOfGenes
        % Initialize Position
        for j=1:numberOfModules
            particle(j).Position=j;
        end
        

        %Initialize Coordinates
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
        rows=1;
        cols=1;
        for j=1:numberOfModules
                particle(j).Irradiation= Irr(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end

        %Current
         rows=1;
        cols=1;
        for j=1:numberOfModules
                particle(j).Current= Imp(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end
        
  
        %Voltages
        rows=1;
        cols=1;
        for j=1:numberOfModules
                particle(j).Voltage= Vmp(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end
        
        %% Make the computation
         
            for k=1:numberOfModules
                a(k)=k;
            end
            idx=randperm(numberOfModules);
            final(i,idx)=a(:);
            
            for j=1:numberOfModules
                tempPosition(j)=final(i,j);
            end
            
            for j=1:numberOfModules
                tempIrradiation(j)=particle(tempPosition(j)).Irradiation;
            end
            
             for j=1:numberOfModules
                tempCurrent(j)=particle(tempPosition(j)).Current;
             end
            
            for j=1:numberOfModules
                tempCoordinates(j)=particle(tempPosition(j)).Coordinates;
            end
            
            %COPY TO PARTICLES
            for j=1:numberOfModules
                particle(j).Position=tempPosition(j);
            end
            gen(i).Number=[particle(:).Position];
            for j=1:numberOfModules
                particle(j).Irradiation=tempIrradiation(j);
            end
            
             for j=1:numberOfModules
                particle(j).Current=tempCurrent(j);
             end
            
            for j=1:numberOfModules
                particle(j).Coordinates=tempCoordinates(j);
            end
            
          %Evaluation
          %Find row columns
          Icolumn(1)= particle(1).Current+particle(2).Current+particle(3).Current;
          Icolumn(2)= particle(4).Current+particle(5).Current+particle(6).Current;
          Icolumn(3)= particle(7).Current+particle(8).Current+particle(9).Current; 
            
          %findmaxRowCurrent
         maxRowCurrent=Icolumn(1);
         for rep=2:3
            if(Icolumn(rep)>maxRowCurrent)
                maxRowCurrent=Icolumn(rep);
            end
         end
         
         %Compute Cost
         gen(i).Cost = CostFunction(Icolumn(1),Icolumn(2),Icolumn(3),maxRowCurrent);
         
         if(gen(i).Cost < BestCost)
             BestCost=gen(i).Cost;
             BestGene=gen(i).Number;
             IcolumnBest(1)= particle(1).Current+particle(2).Current+particle(3).Current;
          IcolumnBest(2)= particle(4).Current+particle(5).Current+particle(6).Current;
          IcolumnBest(3)= particle(7).Current+particle(8).Current+particle(9).Current; 
          
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
         
         %Hold the 20 first Genes
         for i=1:20
             newTempGen(i,:)=tempGen(i,:);
         end
       
         %Create the Mutations
         rowNumber=21;
         for i=1:20
             if(i<=20)
                 
                 randNum1=ceil(10*rand()); %First Index Number
                 while(randNum1==10)
                 randNum1=ceil(10*rand());
                 end;
                 
                 randNum2=ceil(10*rand()); %Second Index Number
                 while(randNum2==10)
                 randNum2=ceil(10*rand());
                 end;
                 
                 %Swap values
                 temporalGen=newTempGen(i,:);
                 tempNumber=temporalGen(randNum1);
                 temporalGen(randNum1)=temporalGen(randNum2);
                 temporalGen(randNum2)=tempNumber;
                 
                 newTempGen(rowNumber,:)=temporalGen;
                 rowNumber=rowNumber+1;
%              else if(i>10 && i<=20)
%                      
%                  randNum1=ceil(10*rand()); %First Index Number
%                  while(randNum1==10)
%                  randNum1=ceil(10*rand());
%                  end;
%                  randNum11=randNum1+1;
%                  if(randNum11==10)
%                      randNum11=randNum1-1;
%                  end
%                  
%                  randNum2=ceil(10*rand()); %Second Index Number
%                  while(randNum2==10)
%                  randNum2=ceil(10*rand());
%                  end;
%                  randNum22=randNum2+1;
%                  while(randNum22>=10 ||randNum22==randNum11)
%                      if(randNum22>=10)
%                       randNum22=randNum2-1;
%                      end
%                      if(randNum22==randNum11)
%                        randNum22=ceil(10*rand());
%                      end
%                  end
%                  
%                  
%                  
%                  temporalGenDouble=newTempGen(i,:);
%                  tempNumberDouble1=temporalGenDouble(randNum1);
%                  tempNumberDouble2=temporalGenDouble(randNum11);
%                  temporalGenDouble(randNum1)=temporalGenDouble(randNum2);
%                  temporalGenDouble(randNum11)=temporalGenDouble(randNum22);
%                  temporalGenDouble(randNum2)=tempNumberDouble1;
%                  temporalGenDouble(randNum22)=tempNumberDouble2;
%                  
%                  
%                  
%                  newTempGen(rowNumber,:)=temporalGenDouble;
%                  rowNumber=rowNumber+1;   
%                  end
             end
                 
                 
         end
         
         %% Create random genes
         for i=1:10
         for k=1:numberOfModules
                aNew(k)=k;
            end
            idxNew=randperm(numberOfModules);
            finalNew(i,idxNew)=aNew(:);
            
            
            newTempGen(rowNumber,:)=finalNew(i,:);
 
            rowNumber=rowNumber+1;
         end
         
         %% Evaluation
         
         for i=1:numberOfGenes
         %Initialize Position
        for v=1:numberOfModules
            particle(v).Position=v;
        end
        

        %Initialize Coordinates
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
        rows=1;
        cols=1;
        for v=1:numberOfModules
                particle(v).Irradiation= Irr(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end

        %Current
         rows=1;
        cols=1;
        for v=1:numberOfModules
                particle(v).Current= Imp(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end
        
  
        %Voltages
        rows=1;
        cols=1;
        for v=1:numberOfModules
                particle(v).Voltage= Vmp(rows,cols);
                rows=rows+1;
            if(rows==4)
            cols=cols+1;
            rows=1;
            end
        end
         
         
         
         
            for j=1:numberOfModules
                tempPosition(j)=newTempGen(i,j);
            end
            
            for j=1:numberOfModules
                tempIrradiation(j)=particle(tempPosition(j)).Irradiation;
            end
            
             for j=1:numberOfModules
                tempCurrent(j)=particle(tempPosition(j)).Current;
             end
            
            for j=1:numberOfModules
                tempCoordinates(j)=particle(tempPosition(j)).Coordinates;
            end
            
            %COPY TO PARTICLES
            for j=1:numberOfModules
                particle(j).Position=tempPosition(j);
            end
            gen(i).Number=[particle(:).Position];
            for j=1:numberOfModules
                particle(j).Irradiation=tempIrradiation(j);
            end
            
             for j=1:numberOfModules
                particle(j).Current=tempCurrent(j);
             end
            
            for j=1:numberOfModules
                particle(j).Coordinates=tempCoordinates(j);
            end
            
          %Evaluation
          %Find row columns
          Icolumn(1)= particle(1).Current+particle(2).Current+particle(3).Current;
          Icolumn(2)= particle(4).Current+particle(5).Current+particle(6).Current;
          Icolumn(3)= particle(7).Current+particle(8).Current+particle(9).Current; 
            
          %findmaxRowCurrent
         maxRowCurrent=Icolumn(1);
         for rep=2:3
            if(Icolumn(rep)>maxRowCurrent)
                maxRowCurrent=Icolumn(rep);
            end
         end
         
         %Compute Cost
         gen(i).Cost = CostFunction(Icolumn(1),Icolumn(2),Icolumn(3),maxRowCurrent);
         
         if(gen(i).Cost <= BestCost)
             BestCost=gen(i).Cost;
             BestGene=gen(i).Number;
             IcolumnBest(1)= particle(1).Current+particle(2).Current+particle(3).Current;
          IcolumnBest(2)= particle(4).Current+particle(5).Current+particle(6).Current;
          IcolumnBest(3)= particle(7).Current+particle(8).Current+particle(9).Current; 
          
          for h=1:numberOfModules
              particle(h).Best.Irradiation=particle(h).Irradiation;
          end
             
         end
         end
 

          [value1,index1]=sort([gen(:).Cost]);
      end    
         %% ShowResults
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
    text(Vtot(PeakIdx)+10, Peak-10, sprintf('GA-GP=%6.1f W', Peak))
    hold on
      