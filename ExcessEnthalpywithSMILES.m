clc
clear 
%% Read data 
T = readtable('ExcessEnthalpyData.xlsx', 'Sheet', 'Fed to network', 'ReadVariableNames',true); % change sheet to include certain functional groups as the main site of data collection 
%  SMILES vectors
[smiles, rowssmiles, colssmiles, lengthssmiles] = SMILES_Encode('SMILES.xlsx',1);
FunctionalGroup1=zeros(length(T.FunctionalGroup1),rowssmiles); %(row, column, page)%, colssmiles,, for matrix
FunctionalGroup2=zeros(length(T.FunctionalGroup1),rowssmiles); % , colssmiles for matrix 

for i=1:length(T.FunctionalGroup1)
    %functional group 1
    key1 = char(cell2mat([T.FunctionalGroup1(i) num2str(T.Chainlength1(i))]));
    key2 = char(cell2mat([T.FunctionalGroup2(i) num2str(T.Chainlength2(i))]));
    %concatenate name and chainlength from data 
    FunctionalGroup1(i,:)=smiles(key1); 
    FunctionalGroup2(i,:)=smiles(key2);
end
%turn these into columns
%FunctionalGroup1(
% C = permute(A,[1 3 2]);
% C = reshape(C,[],size(A,2),1)
 
Tarray = table2array(T(:,[3:11,14])); 
Tarray(:,11)=(Tarray(:,6)./Tarray(:,5)); %temp/ pressure 
Tarray(:,[12,13])=table2array(T(:,[12,13]));%standard enthalpies
Tarray(:,14)= Tarray(:,10)./Tarray(:,6).*Tarray(:,9);%standard enthalpy/temp,
Tdata = [FunctionalGroup1 FunctionalGroup2 Tarray(:,[1:9,12:14])];%array
%for the NN app
m = size(Tdata,2);
% change this to change functional groups included
%inputs to 2D matrix
% A = Tdata(:,1:m,1:6); 
% C = permute(A,[1,3,2]);
% D = reshape(C,length(T.FunctionalGroup1), []);
%inputs and outputs
inputs = Tdata;
target = T.Excessenthalpy; 

%% ANN
nodes = 14; % change number of nodes here
rng = 42;
net = feedforwardnet([nodes,2], 'trainbr'); % type of ANN and regularization  
%choose division of data 
net.divideFcn = 'dividerand'; % interleaved indices used
net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio = 0;
net.divideParam.testRatio = 0.3;
 %  preprocessing as a std dev of 1 and mean of 0 
net.input.processFcns = {'mapminmax'};
net.output.processFcns = {'mapminmax'};
%choose performance function 
net.performFcn = 'mse'; %mean squared error
%choose number of epochs, default = 1000
net.trainParam.epochs = 2000;
% plot fcns 
net.plotFcns = {'plotperform','plottrainstate','ploterrhist','plotregression', 'plotfit','plotwb'};
%choose transfer function in the hidden layer 
net.layers{1}.transferFcn = 'logsig';
%train ANN

[net,tr] = train(net, inputs.', target.'); 
%predictions
y = net(inputs.'); 
error = gsubtract(target.' ,y);
%Performance 
mse=perform(net,target.',y)
trainTarget = (target.').*tr.trainMask{1};
testTarget = (target.').*tr.testMask{1};
valTarget = (target.').*tr.valMask{1};
trainmse = perform(net, trainTarget, y)
testmse = perform(net,testTarget,y)
msediffperc = abs(trainmse-testmse)/trainmse*100
Cyt = corrcoef(y,target.');
R=sqrt(Cyt(2,1));

%% results 
netResults = net(inputs.');
training = netResults(tr.trainInd);
validation = netResults(tr.valInd);
test = netResults(tr.testInd);
in = inputs.';
trainTarget = in(tr.trainInd);
validationTarget = in(tr.valInd);
testTarget = in(tr.testInd);
%errors
errorLim = 5; %setting error limit 
errorlocation = (abs(error)>errorLim)'; %location of errors
errorSize = (error(abs(error)>errorLim))';%magnitude of errors found 
%prediction
%Ypred = predict(net, Tdata(1:12,1:13));

%% RI
wb = formwb(net,net.b,net.iw,net.lw);
[b,IW,LW] = separatewb(net,wb);
r = zeros(nodes,m);
c = zeros(nodes,m);
S = zeros(m,1);
RI = zeros(m,1);
for i=1:m
    for j = 1:nodes
        c(j,i) = IW{1}(j,i)*LW{2,1}(j); 
    end  
end 
for i = 1:m
    for j = 1:nodes
        r(j,i)= abs(c(j,i))/sum(abs(c(j,:)));
    end 
    S(i)=sum(r(:,i));
end
for i = 1:m
    RI(i) = S(i)/sum(S)*100;
end 
%% Prediction plot
clf
indices = 1:20; %indices of the data being plotted 
indicesreal = 1:20;
x = inputs.';
            %plots uncomment for one check of a neural net

indicesopt = {[592:627], [3079:3093],[3156:3175],[3137: 3155]};
%ethane & butan-1-ol at 348.15K, 5000kPa
%hexane and butan-2-ol at 101kpa and 298.15K
%heptane and dodecane at 110kPa and 298.15K
% heptan-1-ol and octan-1-ol at 298.15K and 101kPa
n=2;
tiledlayout(n,size(indicesopt,2)/n)
for i=1:size(indicesopt,2)
    nexttile
    indices= indicesopt{i};
    y1 =T.Excessenthalpy(indices)+T.idealenthalpy(indices);
    y2= y(indices).'+T.idealenthalpy(indices);

    plot(T.Compositioncomponent1(indices), y1,'or',T.Compositioncomponent1(indices), y2, 'k')
    xlabel("Composition component 1 (mol/mol)")
    ylabel("Real Enthalpy (J/mol)")
%     legend("Literature values", "Model")
end
%% Interpolation 
tiledlayout(1,2)
TI = readtable('ExcessEnthalpyData.xlsx', 'Sheet', 'Interpolation');
% One-hot encoding components 
% change length of this vector for more functional groups
[smiles, rowssmiles, colssmiles, lengthssmiles] = SMILES_Encode('SMILES.xlsx',1);
FunctionalGroup1=zeros(length(TI.FunctionalGroup1),rowssmiles); %(row, column, page)%, colssmiles,, for matrix
FunctionalGroup2=zeros(length(TI.FunctionalGroup1),rowssmiles); % , colssmiles for matrix 

for i=1:length(TI.FunctionalGroup1)
    %functional group 1
    key1 = char(cell2mat([TI.FunctionalGroup1(i) num2str(TI.Chainlength1(i))]));
    key2 = char(cell2mat([TI.FunctionalGroup2(i) num2str(TI.Chainlength2(i))]));
    %concatenate name and chainlength from data 
    FunctionalGroup1(i,:)=smiles(key1); 
    FunctionalGroup2(i,:)=smiles(key2);
end
%turn these into columns
%FunctionalGroup1(
% C = permute(A,[1 3 2]);
% C = reshape(C,[],size(A,2),1)
 
TarrayI = table2array(TI(:,[3:11,14])); 
TarrayI(:,11)=(TarrayI(:,6)./TarrayI(:,5)); %temp/ pressure 
TarrayI(:,[12,13])=table2array(TI(:,[12,13]));%standard enthalpies
TarrayI(:,14)= TarrayI(:,10)./TarrayI(:,6).*TarrayI(:,9);%standard enthalpy/temp,
TdataI = [FunctionalGroup1 FunctionalGroup2 TarrayI(:,[1:9,12:14])];%array

inputsI = TdataI;
yI = sim(net,inputsI.'); 
nexttile
plot(TI.Compositioncomponent1(1:15), TI.Excessenthalpy(1:15),'or', TI.Compositioncomponent1(1:15), yI(1:15), 'k')
xlabel("Composition 1(mol/mol)")
ylabel("Excess Enthalpy (J/mol)")
nexttile
plot(TI.Compositioncomponent1(16:45), TI.Excessenthalpy(16:45),'or', TI.Compositioncomponent1(16:45), yI(16:45), 'k')
xlabel("Composition 1(mol/mol)")
ylabel("Excess Enthalpy (J/mol)")
%% Extrapolation 
clf
TE =readtable('ExcessEnthalpyData.xlsx', 'Sheet', 'Extrapolation');
[smiles, rowssmiles, colssmiles, lengthssmiles] = SMILES_Encode('SMILES.xlsx',1);
FunctionalGroup1=zeros(length(TE.FunctionalGroup1),rowssmiles); %(row, column, page)%, colssmiles,, for matrix
FunctionalGroup2=zeros(length(TE.FunctionalGroup1),rowssmiles); % , colssmiles for matrix 

for i=1:length(TE.FunctionalGroup1)
    %functional group 1
    key1 = char(cell2mat([TE.FunctionalGroup1(i) num2str(TE.Chainlength1(i))]));
    key2 = char(cell2mat([TE.FunctionalGroup2(i) num2str(TE.Chainlength2(i))]));
    %concatenate name and chainlength from data 
    FunctionalGroup1(i,:)=smiles(key1); 
    FunctionalGroup2(i,:)=smiles(key2);
end



TarrayE = table2array(TE(:,[3:11,14])); 
TarrayE(:,11)=(TarrayE(:,6)./TarrayE(:,5)); %temp/ pressure = best

TarrayE(:,[12,13])=table2array(TE(:,[12,13]));
TarrayE(:,14)= TarrayE(:,10)./TarrayE(:,6).*TarrayE(:,9);%standard enthalpy/temp, 

TdataE = [FunctionalGroup1 FunctionalGroup2 TarrayE(:,[1:9,12:14])]; % all features besides ideal enthalpy fed to network
inputsE = TdataE;
yE = sim(net,inputsE.');
%yE = yE1+ TE.Idealenthalpy;
tiledlayout(1,2)
nexttile 
index = 1:13;
plot(TE.Compositioncomponent1(index), TE.Excessenthalpy(index),'or', TE.Compositioncomponent1(index), yE(index), 'k')
xlabel("Composition 1(mol/mol)")
ylabel("Excess Enthalpy (J/mol)")
nexttile 
index = 14:24;
plot(TE.Compositioncomponent1(index), TE.Excessenthalpy(index),'or', TE.Compositioncomponent1(index), yE(index), 'k')
xlabel("Composition 1(mol/mol)")
ylabel("Excess Enthalpy (J/mol)")

%%
function [SMILES_encoded, rows, cols, weight] = SMILES_Encode_matrix(filename)
    %encode smiles vectors from an excel file containing smiles as strings 
    % This function codes binary matrices
    T = readtable(filename);
    SMILES_CHARS = {' ', '(', ')', '1', '=',  'C',  'O'};
    %create map of characters and their indices
    M = containers.Map( SMILES_CHARS, 1:length(SMILES_CHARS));
    SMILES_encoded = containers.Map('KeyType', 'char','ValueType', 'any');

    rows = max(strlength(T.SMILES));
    cols = length(SMILES_CHARS);
    X = zeros(rows,cols); %length(T.SMILES), 
    pad_char = ' ';
    weight = zeros(rows,1);
    for i=1:length(T.SMILES)
        
        string_new = T.SMILES(i);
        weight(i,1) = length(T.SMILES(i));
        % all strings to the same length
        for k = 1:rows-length(T.SMILES(i))
            string_new = strcat(string_new, pad_char); 
        end
        string_new=char(string_new);
        %encode the vector
        X = zeros(rows);
        for j = 1:length(string_new)
            c = string_new(j);
            X(j, M(char(c))) = 1;  %
        end
        %create matlab dictionary to store all encoded SMILES vectors 
        key=char(T.namechainlength(i));
        SMILES_encoded(key) = X; %concatenates name and chainlength for the map 
    end 
end 
function [SMILES_encoded, rows, cols, weight] = SMILES_Encode(filename, ascii)
    %encode smiles vectors from an excel file containing smiles as strings 
    % This function codes column vectors with numbers
    T = readtable(filename);
    SMILES_CHARS = {' ', '(', ')', '1', '=',  'C',  'O'};
     % assign acsii characters to characters used in SMILES vectors
    %create map of characters and their indices
    if ascii ==1
        ascii_SMILES = double(cell2mat(SMILES_CHARS));
        M = containers.Map( SMILES_CHARS, ascii_SMILES);
    else 
        M = containers.Map( SMILES_CHARS, 1:length(SMILES_CHARS));
    end 
    
    SMILES_encoded = containers.Map('KeyType', 'char','ValueType', 'any');
    rows = max(strlength(T.SMILES));
    cols = length(SMILES_CHARS);
    % columns = 1 for no binary marix 
    X = zeros(rows,1); %length(T.SMILES), 
    pad_char = ' ';
    weight = zeros(rows,1);
    for i=1:length(T.SMILES)
        
        string_new = T.SMILES(i);
        weight(i,1) = length(T.SMILES(i));
        
        % all strings to the same length, padded with zeros 
        for k = 1:rows-length(T.SMILES(i))
            string_new = strcat(string_new, pad_char); 
        end
        string_new=char(string_new);
        
        %encode the vector
        
        X = zeros(rows,1);
        for j = 1:length(string_new)
            c = string_new(j);
            X(j,1) = M(char(c));  %
        end
        %create matlab dictionary to store all encoded SMILES vectors 
        key=char(T.namechainlength(i));
        SMILES_encoded(key) = X; %concatenates name and chainlength for the map 
    end 
end
