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
Irr(1,2)=900;
Irr(1,3)=900;

Irr(2,1)=900;
Irr(2,2)=200;
Irr(2,3)=900;

Irr(3,1)=900;
Irr(3,2)=900;
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
         
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DIFFERENTIAL EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%       
         for i=1:numberOfGenes
             xGen=tempGen(i,:);
             yGen=[];
             randomSelNumber1=ceil(rand()*50);
             randomSelNumber2=ceil(rand()*50);
             randomSelNumber3=ceil(rand()*50);
             
             aGen=tempGen(randomSelNumber1,:);
             bGen=tempGen(randomSelNumber2,:);
             cGen=tempGen(randomSelNumber3,:);
             
             for j=1:numberOfModules
                CR=rand()
                F=rand()*0.5
                flagCheck=false
                %Check for crossover propability
                if(CR>0.5)
                    yGen(j)=round(aGen(j)+abs(F*(bGen(j)-cGen(j))))
                    %Check for same values
                    for k=1:j
                        if(yGen(j)==yGen(k) && j~=k)
                            flagCheck=true
                            break;
                        end
                    end
                    %Check for out of range and act if something is wrong
                    flagIncrement=true;
                    flagInternal=false;
                    while(yGen(j)>=numberOfModules+1 || flagCheck==true)
                         flagInternal=false;
                        if(yGen(j)>=numberOfModules)
                            flagIncrement=false
                        end
                        if(yGen(j)<numberOfModules && flagIncrement==true)
                            yGen(j)=yGen(j)+1
                        end
                        
                        if(yGen(j)>1 && flagIncrement==false)
                            yGen(j)=yGen(j)-1
                        end
                        
                        for k=1:j
                                if(yGen(j)==yGen(k) && j~=k)
                                    flagInternal=true
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
                            flagCheck=true
                            break;
                        end
                    end
                    flagIncrement=true;
                    flagInternal=false;
                    while(yGen(j)>=numberOfModules+1 || flagCheck==true)
                         flagInternal=false;
                        if(yGen(j)>=numberOfModules)
                            flagIncrement=false
                        end
                        if(yGen(j)<numberOfModules && flagIncrement==true)
                            yGen(j)=yGen(j)+1
                        end
                        
                        if(yGen(j)>1 && flagIncrement==false)
                            yGen(j)=yGen(j)-1
                        end
                        
                        for k=1:j
                                if(yGen(j)==yGen(k) && j~=k)
                                    flagInternal=true
                                    flagCheck=true;
                                end
                        end
                        
                        if(flagInternal==false)
                            flagCheck=false;
                        end     
                    end  
                end   
             end
             %Fill the vectors
             middleTempGen(i,:)=yGen;
         end
  
%          
%          for u=1:50
%              g(u,:)=unique(middleTempGen(u,:))
%          end
%          
%          pop=5;
         
         %Hold the 20 first Genes
         for i=1:20
             newTempGen(i,:)=tempGen(i,:);
         end
       
         %Create the Mutations
         for i=21:40
             newTempGen(i,:)=middleTempGen(i,:);
             
                 
                 
         end
         
         %% Create random genes
         rowNumber=41;
         for i=1:10
         for k=1:numberOfModules
                aNew(k)=k;
            end
            idxNew=randperm(numberOfModules);
            finalNew(i,idxNew)=aNew(:);
            
            
            newTempGen(rowNumber,:)=finalNew(i,:);
 
            rowNumber=rowNumber+1;
         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DIFFERENTIAL EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%           
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
            tempNumber=[particle(:).Position];
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
         tempCost = CostFunction(Icolumn(1),Icolumn(2),Icolumn(3),maxRowCurrent);
         
         if(tempCost < gen(i).Cost)
            tempGen(i,:)=yGen;
            gen(i).Cost=tempCost;
            gen(i).Number=tempNumber;
         end
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
    text(Vtot(PeakIdx), Peak-20, sprintf('DE-GP=%6.1f W', Peak))
    hold off
    leg=legend('PSO','PSO-GP','GA','GA-GP','DE','DE-GP');
    leg.Location=('northwest');
    close_system('tepeRarx.slx');
      