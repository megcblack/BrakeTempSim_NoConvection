%RIT Formula SAE Brake Temperature Model

%Original authors: Steven Reuter, Jon Washington

%Modified 9/10/24 by Ethan Onslow

%Modified 9/24/2025 by Oliver Owen (omo9555@rit.edu,
%oliver.owen.pers@gmail.com), edits:
    % updated data input to use "readcell" and "readmatrix" to improve compatability and reduce runtime.
    % changed variable names to improve readability.
    % updated loop structure to allow for the reading and manipulation of multiple sets of data per iteration.

%Modified 9/29/25 by Meghan Black

%Purpose:
%The purpose of this program is to obtain the temperature of the brake rotor based on vehicle parameters and thermal properties

%% INPUT COLLECTION

clc, clear, format compact %clear variables, command window, and format output
fprintf('Welcome to the RIT Brake Temperature Simulation\n\n')

% READ PARAMS 
Car_Parameter = getCarParams;

% Get Drive Data File Path
driveDataFile = getDriveData;

% Other input
disp('Would you like to generate a pdf?\n')
pdf = input('yes = 0 no = 1');
if pdf == 0
    [pdfFile, pdfPath] = uiputfile('*.pdf', 'Save Data As', "Results_NoConvection.pdf");
    filePath = strcat(pdfPath, pdfFile);
end
disp('would you like suppress matlab gui graph output?')
hide_graphs = input('yes = 0, no = 1');

%Definition of variables 
VehicleMass = Car_Parameter(1); % [kg]
RotorMass = Car_Parameter(2); % [kg]
RotorArea = Car_Parameter(3); % [m^2]
dragCoeff = Car_Parameter(4); % [-]
frontalArea = Car_Parameter(5); % [m^2]
brakebias = Car_Parameter(6); % [-], portion into front axle
rotInertia = Car_Parameter(7); % [kg*m^2]
wheelRad = Car_Parameter(8); % [m]
airDensity = Car_Parameter(9); % [kg/m^3] 

%Import data from Excel
[~, sheetNames] = xlsfinfo(driveDataFile); % create array of sheet names 
numDataSets = numel(sheetNames); % returns length of sheetNames array, i.e., number of sheets/data sets

% initialize and fill array of vectors for each parameter of interest
timeArray = cell(numDataSets, 1 );
brakeTempArrayC = cell(numDataSets, 1);
speedArrayMPH = cell(numDataSets, 1);
brakePressArray = cell(numDataSets, 1);

for n = 1:1:numDataSets
    timeArray{n} = readmatrix(driveDataFile, 'Sheet', n, 'Range', 'A:A'); % [sec]
    brakeTempArrayC{n} = readmatrix(driveDataFile, 'Sheet', n, 'Range', 'B:B'); % [deg C]
    speedArrayMPH{n} = readmatrix(driveDataFile, 'Sheet', n, 'Range', 'C:C'); % [MPH]
    brakePressArray{n} = readmatrix(driveDataFile, 'Sheet', n, 'Range', 'D:D'); % [psi]
end

%% Filtering and unit conversion

%initialize and fill new arrays of vectors for filtered/converted PoI
brakeTempArrayF = cell(numDataSets, 1);
speedArray = cell(numDataSets, 1);

for m = 1:1:numDataSets
    % strip data of headers; assumes the first row is headers
    timeArray{m}(1) = [];
    brakeTempArrayC{m}(1) = []; 
    speedArrayMPH{m}(1) = [];
    brakePressArray{m}(1) = [];

    % correct noise in speed data and convert from MPH to m/s
    for mm = 1:1:length(speedArrayMPH{m})
        if speedArrayMPH{m}(mm) > 80
        speedArrayMPH{m}(mm) = 0;
        end
        speedArray{m}(mm) = speedArrayMPH{m}(mm)*.277778;
    end
    % convert brake temp to Farenheit
    for mn = 1:1:length(brakeTempArrayC{m})
        brakeTempArrayF{m}(mn) = (brakeTempArrayC{m}(mn) * 9/5) + 32;
    end
end

  %% Tuning Parameters
    x1=1.85; %convection coeffecient trend line
    b1=60;
    x2=-.000172; %percent temperature into pads
    b2= 0.595;

for curDataSet = 1:1:numDataSets
    %% Initial Conditions
    TambC{curDataSet} = brakeTempArrayC{curDataSet}(1); % [degC] ambient air temperature, second value in column as data is assumed to have headers
    TambK{curDataSet} = TambC{curDataSet} + 273; % [degK]
    RotorTempArrayK{curDataSet}(1) = TambK{curDataSet}; % [degK]

    %% Calculate displacement
    d{curDataSet}(1) = 0; % [m] intial displacement
    for n = 2:1:length(timeArray{curDataSet}) % 
        timestep = timeArray{curDataSet}(n)-timeArray{curDataSet}(n-1); % [sec] time step to increase the D
        avgspd = (speedArray{curDataSet}(n-1) + speedArray{curDataSet}(n))/2 ; % [m/s] average speed over the timestep
        d{curDataSet}(n) = avgspd*timestep; % [m]
    end

    %% Calculate air braking energy 
    for k = 2:1:length(timeArray{curDataSet})
        Fdrag{curDataSet}(k) = 0.5 * airDensity * (speedArray{curDataSet}(k))^2 * dragCoeff * frontalArea; % [N] instantaneous drag force
        Edrag{curDataSet}(k) = (Fdrag{curDataSet}(k) * d{curDataSet}(k)); % [N*m] instantaneous air drag energy
    end

    %% Calculate rotor temp
    for i = 2:1:length(timeArray{curDataSet}) 
        prevSpeed = speedArray{curDataSet}(i-1); %Previous Linear Speed m/s
        newSpeed = speedArray{curDataSet}(i); %New Linear Speed m/s
        omegaP = prevSpeed/wheelRad; %Previous Rotational Speed rad/s
        omegaN = newSpeed/wheelRad; %New Rotational Speed rad/s
        DS = (newSpeed - prevSpeed); %Speed delta m/s
        prevTemp{curDataSet}(i) = RotorTempArrayK{curDataSet}(i-1); %stores new initial temp
        
        tbrake = timeArray{curDataSet}(i)-timeArray{curDataSet}(i-1); %time to current loop step
        h_w=x1*speedArray{curDataSet}(i) + b1;
        h{curDataSet}(i)=h_w/1000;%kW/(m^2 K)
        Rrotor=1/(h{curDataSet}(i)*RotorArea); %thermal resistance of rotor convection
        RotorDiameter = .206; %m
        PadFrac{curDataSet}(i)=prevTemp{curDataSet}(i)*x2 + b2; %percent braking power into pads
        SpecHeat = 0.0005*prevTemp{curDataSet}(i)+0.2813; %specific heat capacity of 4130 [J/g*K]
        
        if DS < -5 %if speed sensor is tripped, the simulation will not see braking by accident
           DS=1;
        end
    
        if DS < 0 & brakePressArray{curDataSet}(i)> min(brakePressArray{curDataSet})% determines whether you are braking or accelerating, needs brake pressure input
           Energy1{curDataSet}(i) = (.5*VehicleMass*(prevSpeed^2 - newSpeed^2)); %Linear Kinetic Energy J
           Energy2{curDataSet}(i) = 4*(.5*rotInertia*(omegaP^2-omegaN^2)); %Rotational Kinetic Energy J
           Energy{curDataSet}(i) = Energy1{curDataSet}(i) + Energy2{curDataSet}(i); %Total Energy J
           %total energy from braking
           AeroFrac{curDataSet}(i) = Edrag{curDataSet}(i)/(Energy{curDataSet}(i)); %areo dynamic reduction in energy, unitless
           BrakeFrac = .82 ; %based on average engine brake data
           
           CorrectedEnergy{curDataSet}(i) = Energy{curDataSet}(i) * 0.5*brakebias*(1-AeroFrac{curDataSet}(i))*(1-PadFrac{curDataSet}(i))*(BrakeFrac)*0.8; %J --> account for brake bias
           deltaTK = ((CorrectedEnergy{curDataSet}(i)/1000)/(RotorMass * SpecHeat)); %rise in temp in rotor kelvin, energy converted to kJ
           RotorTempArrayK{curDataSet}(i) = deltaTK + prevTemp{curDataSet}(i); %adds change in T to previous temp, divide by two because 
           Power{curDataSet}(i) = CorrectedEnergy{curDataSet}(i)/tbrake;
    
           qout{curDataSet}(i) = (RotorTempArrayK{curDataSet}(i)-TambK{curDataSet})/Rrotor; %Kwatts
           Eout{curDataSet}(i) = qout{curDataSet}(i)*tbrake; %Kjoules
           deltaTKout = ((Eout{curDataSet}(i))/(RotorMass * SpecHeat));
           RotorTempArrayK{curDataSet}(i) = RotorTempArrayK{curDataSet}(i)-(deltaTKout);
    
        else %brakes not applied, cooling
           qout{curDataSet}(i) = (prevTemp{curDataSet}(i)-TambK{curDataSet})/Rrotor; %Kwatts
           Eout{curDataSet}(i) = qout{curDataSet}(i)*tbrake; %Kjoules
           deltaTKout = ((Eout{curDataSet}(i))/(RotorMass * SpecHeat)); 
           RotorTempArrayK{curDataSet}(i) = prevTemp{curDataSet}(i)-(deltaTKout);
           Energy1{curDataSet}(i) =0;
           Energy2{curDataSet}(i) = 0;
           Energy{curDataSet}(i) =0;
           CorrectedEnergy{curDataSet}(i) =0;
           EBE{curDataSet}(i) =0;
           Power{curDataSet}(i) = CorrectedEnergy{curDataSet}(i)/tbrake;
        end
    end

    % convert rotor temp from K to F, fill last value of array, and obtain maxium temperature
    for c = 1:1:length(RotorTempArrayK{curDataSet})
        RotorTempArrayF{curDataSet}(c) = RotorTempArrayK{curDataSet}(c)*(9/5)-459.67;
    end
    RotorTempArrayF{curDataSet}(length(timeArray{curDataSet}))= RotorTempArrayF{curDataSet}(length(timeArray{curDataSet}) - 1);
    PeakTemp{curDataSet} = max(RotorTempArrayF{curDataSet});
    fprintf('\nPeak Temperature for Data Set %d: %3.2f degF', curDataSet, PeakTemp{curDataSet})

    %% Percent difference
    PercErr{curDataSet} = zeros(length(RotorTempArrayF),1);
    for n = 1:1:length(RotorTempArrayF{curDataSet})
        PercErr{curDataSet}(n)=(abs(RotorTempArrayF{curDataSet}(n)-brakeTempArrayF{curDataSet}(n))/brakeTempArrayF{curDataSet}(n))*100;
        if PercErr{curDataSet}(n) > 100
            PercErr{curDataSet}(n) = PercErr{curDataSet}(n-1);
        else
            PercErr{curDataSet}(n) = PercErr{curDataSet}(n);
        end
    end
    AvgP{curDataSet} = mean(PercErr{curDataSet});
    fprintf('\nPercent Error of Temperature for Data Set %d = %3.2f%%\n ',curDataSet, AvgP{curDataSet})        

    %% Plots
    f = figure('Name', "Brake Temperature Comparison " + num2str(curDataSet)); % create temperature comparison figure
    title(['Brake Temperature Comparison for Data Set ', num2str(curDataSet)])
    hold on
    plot(timeArray{curDataSet},RotorTempArrayF{curDataSet},'-r') %sim temps, red line.
    plot(timeArray{curDataSet},brakeTempArrayF{curDataSet},'-b') %real temps, blue line
    ylabel('Temperature (Deg F)') 
    legend('Simulated Temperature','Actual Temperature','Location','northwest')
    xlabel('Time (seconds)')
    if pdf == 0
        exportgraphics(f.CurrentAxes, filePath, 'Append', true);
    end
    
    f = figure('Name', "Temperature Percent Error " + num2str(curDataSet));
    plot(timeArray{curDataSet},PercErr{curDataSet})
    title(['Temperature Percent Error for Data Set ', num2str(curDataSet)])
    xlabel('Time (seconds)')
    ylabel('% Error')
    if pdf == 0
        exportgraphics(f.CurrentAxes, filePath, 'Append', true);
    end
    
    f = figure('Name', "Power Over Time " + num2str(curDataSet));
    plot(timeArray{curDataSet},Power{curDataSet})
    title(['Power Over Time for Data Set ', num2str(curDataSet)])
    xlabel('Time (seconds)')
    ylabel('Power (Watts)')
    if pdf == 0
        exportgraphics(f.CurrentAxes, filePath, 'Append', true);
    end
    f = figure('Name', "Pad Fraction Over Time " + num2str(curDataSet));
    plot(timeArray{curDataSet},PadFrac{curDataSet})
    title(['Pad Fraction Over Time for Data Set ', num2str(curDataSet)])
    xlabel('Time (seconds)')
    ylabel('Pad Fraction')
    if pdf == 0
        exportgraphics(f.CurrentAxes, filePath, 'Append', true);
    end

    if hide_graphs == 0
        close all
    end
end



function Car_Parameter = getCarParams()
% GETCARPARAMSFILE - Gets the Car Parameter Data File. Defaults to
% 'CarParameters.xlsx' or asks the user to enter new path
%
% Input Arguements:
%   file path - when prompted, the file path must have the car parameters in
%   column G in the following order:
%   1) Vehicle Mass
%   2) Specific Heat
%   2) Rotor Mass
%   3) Rotor Surface Area
%   4) Drag Coefficient
%   5) Frontal Area
%   6) Brake Bias
%   7) Rotational Inertia of Wheel
%   8) Wheel Radius
%   9) Air Density
%  10) Pedal Ratio
%  11) Master Cylinder Size
%  12) # Front Caliper Pistons
%  13) # Rear Caliper Pistons
%  14) Front Piston Diameter
%  15) Rear Piston Diameter
%
% Output Arguements:
%   Car_Parameter = Car Parameter array 
%
    carParamNames = ["Vehicle Mass", "Specific Heat", "Rotor Mass", "Rotor Surface Area", "Drag Coefficient", "Frontal Area", "Brake Bias", "Rotational Intertia of Wheel", "Wheel Radius", "Air Density", "Pedal Ratio", "Master Cylinder Size", "Front Caliper Pistons", "Rear Caliper Pistons", "Front Piston Diameter", "Rear Piston Diameter"]';
    Units = ["kg", "J/g*C", "kg", "m^2", "-", "m^2", "-","kg*m^2", "m", "kg*m^3", "-", "m", "-", "-", "m", "m"]';
    Units = categorical(Units);
    disp('The following car parameters are being used by the simulation')
    Car_Parameter = readmatrix('CarParameters.xlsx','Sheet', 1,'Range', 'G:G');
    T = table(Car_Parameter,Units,'RowNames', carParamNames);
    disp(T)
    
    disp('would you like to change these values?')
    a = input('yes = 0 no = 1');
    if a == 0
        winopen("CarParameters.xlsx");
        Car_Parameter = getCarParams;
    end
       
end


function driveDataFile = getDriveData()
% DRIVEDATAFILE - read the drive data file
%
% Output Arguements: 
%   driveDataFile - drive data
%
    driveDataFile = readcell('CarParameters.xlsx','Sheet', 1, 'Range', 'J3:J3');
    fprintf('The current data file to be read is the following\n')
    disp(driveDataFile)
    driveDataFile=driveDataFile{1}; % cast driveDataFile from cell type to string type to use as input to readmatrix

    disp('would you like to change these values?')
    a = input('yes = 0 no = 1');
    if a == 0
        winopen("CarParameters.xlsx");
        driveDataFile = getDriveData;
    end
end







