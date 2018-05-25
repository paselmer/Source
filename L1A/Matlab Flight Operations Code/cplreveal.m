function varargout = cplreveal(varargin)
% CPLREVEAL MATLAB code for cplreveal.fig
%      CPLREVEAL, by itself, creates a new CPLREVEAL or raises the existing
%      singleton*.
%
%      H = CPLREVEAL returns the handle to a new CPLREVEAL or the handle to
%      the existing singleton*.
%
%      CPLREVEAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CPLREVEAL.M with the given input arguments.
%
%      CPLREVEAL('Property','Value',...) creates a new CPLREVEAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cplreveal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cplreveal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cplreveal

% Last Modified by GUIDE v2.5 06-Dec-2017 13:41:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cplreveal_OpeningFcn, ...
                   'gui_OutputFcn',  @cplreveal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before cplreveal is made visible.
function cplreveal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cplreveal (see VARARGIN)

% Choose default command line output for cplreveal
handles.output = hObject;

% Update handles structure and Import the data
% [handles.newdumpcpl,handles.statuscpl]=urlread('https://asp-interface.arc.nasa.gov/API/binary_packet_data/N809NA/CPL_BIN');
% [handles.newdumpcats,handles.statuscats]=urlread('http://asp-interface-2.arc.nasa.gov/API/binary_packet_data/N806NA/CATS_BIN');
% [handles.newdumpiwg,handles.statusiwg]=urlread('https://asp-interface.arc.nasa.gov/API/parameter_data/N809NA/IWG1');
 handles.wgetpath = 'C:\Users\akupchoc.NDC\Documents\MobaXterm\slash\bin\wget.exe';
 handles.currentdir = cd;
 handles.iwg1path = horzcat(handles.currentdir,filesep,'iwg1');
 handles.cplpath = horzcat(handles.currentdir,filesep,'cpl');
 handles.camalpath = horzcat(handles.currentdir,filesep,'camal');
 handles.cplurl='https://asp-interface.arc.nasa.gov/API/binary_packet_data/N806NA/CPL_BIN_IN';
 %handles.catsurl='http://asp-interface-2.arc.nasa.gov/API/binary_packet_data/N806NA/CATS_BIN';
 handles.camalurl='https://asp-interface.arc.nasa.gov/API/binary_packet_data/N806NA/CAMAL_INMARSAT';
 handles.iwg1url='https://asp-interface.arc.nasa.gov/API/parameter_data/N806NA/IWG1';
handles.newdumpcpl = zeros(1,1680);
handles.statuscpl = 1;
handles.cplprof=zeros(360,833);
handles.camalprof=zeros(360,1000);
handles.currentdirectory = cd;
handles.cpl=zeros(1,833);
handles.newTimeiwg='';
handles.newTimecamal='';
zero4rainbow = zeros(1,21);
one4rainbow = ones(1,21);
taper4rainbow = (0:.05:1);
handles.UnixStartTime = [1970 1 1 0 0 0];
handles.GPSTime = [1980 1 6 0 0 0];
%rainbow = [zero4rainbow' zero4rainbow' taper4rainbow'; zero4rainbow' taper4rainbow' one4rainbow'; zero4rainbow' one4rainbow' 1-taper4rainbow'; taper4rainbow' one4rainbow' zero4rainbow'; one4rainbow' 1-taper4rainbow' zero4rainbow'; 1 1 1];

rainbow =[         0         0         0;    
    0.0500      0.02745098     0.1000;
    0.1300      0.02745098     0.2000;
    0.1800      0.02745098     0.3200;
    0.2300      0.02745098     0.3600;
    0.2800      0.02745098     0.4000;
    0.3300      0.02745098     0.4400;
    0.3800      0.02745098     0.4900;
    0.3843      0.02745098     0.4941;
    0.3500      0.02745098     0.4941;
    0.3400      0.02745098     0.5000;
    0.3100      0.02745098     0.5300;
    0.2500      0.02745098     0.5700;
    0.2000      0.02745098     0.6200;
    0.1500      0.02745098     0.7000;
    0.1000      0.02745098     0.8500;
    0.0500      0.02745098     0.9400;
         0         0    1.0000;
         0    0.0500    1.0000;
         0    0.1000    1.0000;
         0    0.1500    1.0000;
         0    0.2000    1.0000;
         0    0.2500    1.0000;
         0    0.3000    1.0000;
         0    0.3500    1.0000;
         0    0.4000    1.0000;
         0    0.4500    1.0000;
         0    0.5000    1.0000;
         0    0.5500    1.0000;
         0    0.6000    1.0000;
         0    0.6500    1.0000;
         0    0.7000    1.0000;
         0    0.7500    1.0000;
         0    0.8000    1.0000;
         0    0.8500    1.0000;
         0    0.9000    1.0000;
         0    0.9500    1.0000;
         0    1.0000    1.0000;
         0    1.0000    1.0000;
         0    1.0000    0.9500;
         0    1.0000    0.9000;
         0    1.0000    0.8500;
         0    1.0000    0.8000;
         0    1.0000    0.7500;
         0    1.0000    0.7000;
         0    1.0000    0.6500;
         0    1.0000    0.6000;
         0    1.0000    0.5500;
         0    1.0000    0.5000;
         0    1.0000    0.4500;
         0    1.0000    0.4000;
         0    1.0000    0.3500;
         0    1.0000    0.3000;
         0    1.0000    0.2500;
         0    1.0000    0.2000;
         0    1.0000    0.1500;
         0    1.0000    0.1000;
         0    1.0000    0.0500;
         0    1.0000         0;
         0    1.0000         0;
    0.0500    1.0000         0;
    0.1000    1.0000         0;
    0.1500    1.0000         0;
    0.2000    1.0000         0;
    0.2500    1.0000         0;
    0.3000    1.0000         0;
    0.3500    1.0000         0;
    0.4000    1.0000         0;
    0.4500    1.0000         0;
    0.5000    1.0000         0;
    0.5500    1.0000         0;
    0.6000    1.0000         0;
    0.6500    1.0000         0;
    0.7000    1.0000         0;
    0.7500    1.0000         0;
    0.8000    1.0000         0;
    0.8500    1.0000         0;
    0.9000    1.0000         0;
    0.9500    1.0000         0;
    1.0000    1.0000         0;
    1.0000    1.0000         0;
    1.0000    0.9500         0;
    1.0000    0.9000         0;
    1.0000    0.8500         0;
    1.0000    0.8000         0;
    1.0000    0.7500         0;
    1.0000    0.7000         0;
    1.0000    0.6500         0;
    1.0000    0.6000         0;
    1.0000    0.5500         0;
    1.0000    0.5000         0;
    1.0000    0.4500         0;
    1.0000    0.4000         0;
    1.0000    0.3500         0;
    1.0000    0.3000         0;
    1.0000    0.2500         0;
    1.0000    0.2000         0;
    1.0000    0.1500         0;
    1.0000    0.1000         0;
    1.0000    0.0500         0;
    1.0000         0         0;
    1.0000    1.0000    1.0000];


handles.rainbow = rainbow;
handles.colors={'Rainbow','Jet','HSV','Hot', 'Cool','Spring','Summer','Autumn','Winter'...
             'Gray','Bone','Copper','Pink'};
guidata(hObject, handles);
% This sets up the initial plot - only do when we are invisible
% so window can get raised using cplreveal.
if strcmp(get(hObject,'Visible'),'off')
    membrane(1);
end

% UIWAIT makes cplreveal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cplreveal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
if strcmp(get(handles.pushbutton1,'String'),'Start')
    set(handles.pushbutton1,'String', 'Stop');
    while strcmp(get(handles.pushbutton1,'String'),'Stop')
        CompTime = clock;
        set(handles.text2,'String', sprintf('Computer Time %02i:%02i:%i  %02i:%02i:%02.2i',...
            CompTime(2),CompTime(3),CompTime(1), CompTime(4),CompTime(5),round(CompTime(6))));
        %%
        %CAMAL Update
        
%         valuefromslider = sprintf('The slider is reporting %i ',get(handles.slider1, 'Value'));
%         display(valuefromslider);
        
        
        if strcmp(get(handles.uitable2,'Visible'),'on')
            
            system(horzcat(handles.wgetpath,' wget -q -O ',handles.camalpath,' -r -nd --no-check-certificate ',handles.camalurl));        %not tested with CAMAL 11/23/2017
            
            flun = fopen('camal');
            handles.newdumpcamal  = fread(flun,2287,'uint8');                       %269hk + 18 + 2*bins       
            
            fclose(flun);
            
            if(isempty(handles.newdumpcamal))
                display('There is currently no CAMAL data');
            else
            file= horzcat(handles.currentdirectory,filesep, 'camal');        %%%  
            flun = fopen(file, 'r');
            hkcount=1;
            frewind(flun);
            
            %Record Header
            APID=fread(flun,1,'uint16','b');
            SequenceCount=fread(flun,1,'uint16','b');
            ByteCount=fread(flun,1,'uint16','b');
            UnixTime=fread(flun,1,'uint32');
            SubSeconds=fread(flun,1,'uint16');
            GoodCommandCount=fread(flun,1,'int32');
            BadCommandCount=fread(flun,1,'int32');
            %MCS 24
            MCSStatus=fread(flun,1,'int8');
            MCSState=fread(flun,1,'int32');
            NumMCSPacketsCollected=fread(flun,1,'int32');
            %Control
            ControlState=fread(flun,1,'int32');
            ControlSubmode=fread(flun,1,'int32');
            %Scanner
            ScannerState=fread(flun,1,'int32');
            ScannerPos = fread(flun,1,'float32');
            ScannerTemps(1:16) = fread(flun,16,'float32');
            %A/D and discrete
            IOinputs=fread(flun,1,'int8');
            IOoutputs = fread(flun,1,'int8');
            ADvoltages(1:32) = fread(flun,32,'float32');
            %GPS
            GpsWeek=fread(flun,1,'uint32');
            GpsMilliseconds=fread(flun,1,'uint32');
            Latitude = fread(flun,1,'float32');
            Longitude = fread(flun,1,'float32');
            Altitude = fread(flun,1,'float32');
            NumGPSbytesread = fread(flun,1,'uint32');
            %NASDAT
            NASDATflags=fread(flun,1,'uint8');
            NumIWG1pktsrec=fread(flun,1,'uint32');
            Flags=fread(flun,1,'uint8');
            
            %Shapshot
            MCSTimestamp=fread(flun,1,'uint32');
            MCSTimestampMilli=fread(flun,1,'uint16');            
            ChannelIndex=fread(flun,1,'int32');
            NumberofBins=fread(flun,1,'int32');
            Binwidth_sec=fread(flun,1,'float32');           
            newcamalprof=fread(flun,NumberofBins,'uint16');
            handles.camalprof(2:360,:)=handles.camalprof(1:359,:);
            handles.camalprof(1,:) = newcamalprof;
            
             
            ScannerTemps(2:end) = (ScannerTemps(2:end).^2*(-1.9724214))+(ScannerTemps(2:end).*(120.992655))-443.417057; 
            ScannerTemps(1) = (ScannerTemps(1).*(7.4134))-3.20342;            
            ADvoltages(24) = (ADvoltages(:,24).*(7.4134))-3.20342;
            ADvoltages(25:26) =  (ADvoltages(:,25:26).^2*(-1.9724214))+(ADvoltages(:,25:26).*(120.992655))-443.417057 + 19;          %Added 19 degree offset for box temps
            Longitude = Longitude/100;
            Latitude = Latitude/100;
            UnixTime = datevec(datenum(UnixTime(:)/86400)+datenum(handles.UnixStartTime));
            Time = datevec(datenum(datenum(handles.GPSTime)+7*(GpsWeek(:))+(GpsMilliseconds(:)/86400000)));
            
            fclose(flun);
            delete('camal');
            %camalpacketlength = sprintf('CAMAL packet length is %i bytes long',length(handles.newdumpcamal));
            %display(camalpacketlength);

              CAMAL=num2cell([APID; SequenceCount; ByteCount; datenum(UnixTime); SubSeconds; GoodCommandCount; BadCommandCount; MCSStatus; MCSState; NumMCSPacketsCollected; ControlState; ControlSubmode; ...
                ScannerState; ScannerPos; ScannerTemps(1); ScannerTemps(2); ScannerTemps(3); ScannerTemps(4); ScannerTemps(5); ScannerTemps(6); ScannerTemps(7); ScannerTemps(8); ScannerTemps(9); ...
                ScannerTemps(10); ScannerTemps(11); ScannerTemps(12); ScannerTemps(13); ScannerTemps(14); ScannerTemps(15); ScannerTemps(16);  IOinputs; IOoutputs; ...
                ADvoltages(1); ADvoltages(2); ADvoltages(3); ADvoltages(4); ADvoltages(5); ADvoltages(6); ADvoltages(7); ADvoltages(8); ADvoltages(9); ADvoltages(10); ADvoltages(11); ...
                ADvoltages(12); ADvoltages(13); ADvoltages(14); ADvoltages(15); ADvoltages(16); ADvoltages(17); ADvoltages(18); ADvoltages(19); ADvoltages(20); ADvoltages(21); ...
                ADvoltages(22); ADvoltages(23); ADvoltages(24); ADvoltages(25); ADvoltages(26); ADvoltages(27); ADvoltages(28); ADvoltages(29); ADvoltages(30); ADvoltages(31); ADvoltages(32); ...
                GpsWeek; GpsMilliseconds; Latitude; Longitude; Altitude; NumGPSbytesread; NASDATflags; NumIWG1pktsrec; Flags; datenum(Time)]);
            
            switch MCSState
                case 1
                    CAMAL{9} = 'NOT_OPENED';
                case 2
                    CAMAL{9} = 'IDLE';
%                 case 3
%                     CAMAL{9} = 'START_COLLECTING';
                case 3
                    CAMAL{9} = 'COLLECTING';
                case 4
                    CAMAL{9} = 'STOP_COLLECTING';
            end;            
            
            switch ControlState
           %     case 1
           %         CAMAL{11} = 'INITIALIZING';
                case 1
                    CAMAL{11} = 'IDLE';
                case 2
                    CAMAL{11} = 'START_SYSTEM';
                case 3
                    CAMAL{11} = 'WAIT_FOR_LASER_ENABLE';
                case 4
                    CAMAL{11} = 'MOVE_TO_NEXT_SCANNER_POSITION';
                case 5
                    CAMAL{11} = 'COLLECT_DATA_AT_POSITION';
                case 6    
                    CAMAL{11} = 'COLLECT_DATA';

            end;
            switch IOinputs
                case 2
                    CAMAL{31} = 'Sense Laser Enable';
                case 4
                    CAMAL{31} = 'Sense Scanner Enable';
                case 6
                    CAMAL{31} = 'Sense Laser and Scanner Enable';

            end;
            switch ScannerState
                case 0
                    CAMAL{13} = 'SCANNER_NOT_FOUND';
                case 1
                    CAMAL{13} = 'IDLE';
                case 2
                    CAMAL{13} = 'MOVING';
            end;
            switch Flags
                case 0
                    CAMAL{73} = 'Data Drive NOT Seen';
                case 1
                    CAMAL{73} = 'DATA DRIVE SEEN';
            end;
                    
            set(handles.uitable2, 'Data',CAMAL);
            if(strcmp(handles.newTimecamal,num2str(GpsMilliseconds)))
                set(handles.uitable2,'ForegroundColor',[0.31, 0.31, 0.31]);
            else
                set(handles.uitable2,'ForegroundColor',[0, 0, 0]);
            end;
            handles.newTimecamal=num2str(GpsMilliseconds);
            if (get(handles.popupmenu1, 'Value')>= 6)           %was 7
                axes(handles.axes2);
                if get(handles.radiobutton2,'Value')
                    rotate3d on;
                    pause(5);
                else
                    rotate3d off;
                end;    
                    clr=get(handles.popupmenu4,'Value');
                    if (clr > 1)
                        colormap(lower(handles.colors{clr}));
                    else
                        colormap(handles.rainbow)
                    end;
                    % set(hObject, 'String', {'Image', 'Surf','Mesh','Mesh Contour (meshc)', 'Mesh(meshz)','Waterfall',...
                    %     'Contour Fill (contourf)','Contour 3D (contour3)', 'Surf Contour (surfc)','Ribbon' });
                    switch get(handles.popupmenu3, 'Value');
                        case 1
                            maxuse = (get(handles.slider1, 'Value')*1000)+10;
                            imagesc(handles.camalprof',[0 maxuse]);
                        case 2
                            surf(handles.camalprof);
                        case 3
                            mesh(handles.camalprof);
                        case 4
                            meshc(handles.camalprof);
                        case 5
                            meshz(handles.camalprof);
                        case 6
                            contourf(handles.camalprof);
                        case 7
                            contour3(handles.camalprof);
                        case 8
                            surfc(handles.camalprof);
                        case 9
                            ribbon(handles.camalprof);
                        case 10
                            pcolor(handles.camalprof);
                    end;
                    if get(handles.radiobutton1,'Value')
                        axes(handles.axes2);
                        if strcmp(get(colorbar,'Visible'),'off');
                            set(colorbar,'Visible','on');
                        end;
                    else
                        axes(handles.axes2);
                        colorbar('hide');
                    end
                axes(handles.axes1);
                
                %%
            axes(handles.axes1);
            popup_sel_index = get(handles.popupmenu1, 'Value');
            if (popup_sel_index >=5)
                set(handles.edit1, 'String', GPSAltitude);
            end;
            handles.camalprof(1,:) = handles.camalprof(1,:) - mean(handles.camalprof(1,800:1000));             %ADD ANOTHER RADIO BUTTON TO TOGGLE BACKGROUND SUB
            switch popup_sel_index      %'Raw', 'Background Subtracted','Raw - Altitude', 'Background Subtracted - Altitude'
                case 1
                c = plot(newcamalprof); xlim([0 833]); ylim([0 400]);
                set(c,'YData',newcamalprof);
                    pause(1);
                case 2
                    newcamalprof(1,:) = newcamalprof(1,:) - mean(newcamalprof(1,800:1000)); 
                    plot(newcamalprof); xlim([0 833]); ylim([0 400]);
                    pause(1);
                case 3
                    if strcmp('NaN',num2str(str2double(get(handles.edit1,'String'))));
                        display('The Altitude must be numerical')
                        pause(1);
                    else
                        plot(newcamalprof); xlim([0 1000]); ylim([0 400]);
                        set(gca,'XTickLabel',{str2double(get(handles.edit1,'String'))-((0:100:1000)*30)});
                        pause(1);
                    end;
                case 4
                    if strcmp('NaN',num2str(str2double(get(handles.edit1,'String'))));
                        display('The Altitude must be numerical')
                        pause(1);
                    else
                        beep
                        pause(1);
                    end;
                case 5
                        plot(newcamalprof); xlim([0 833]); ylim([0 200]);
                        set(gca,'XTickLabel',{str2double(get(handles.edit1,'String'))-((0:100:1000)*30)});
                        pause(1);
                case 6
                        newcamalprof = newcamalprof - mean(newcamalprof(800:1000));         %UPDATED
                        plot(newcamalprof); xlim([0 1000]); ylim([0 200]);
                        set(gca,'XTickLabel',{str2double(get(handles.edit1,'String'))-((0:100:1000)*30)});
                        pause(1);
            end;
                %%
%             	plot(newcamalprof);
%                 ylim([1 300]);
%                 xlabel('Bins');
%                 ylabel('Counts');
            end;
            end;
        end;
    %%
    %IWG1 Update
if strcmp(get(handles.uitable3,'Visible'),'on')
    %[handles.newdumpiwg,handles.statusiwg]=urlread(handles.iwg1url);
    system(horzcat(handles.wgetpath,' wget -q -O ',handles.iwg1path,' -r -nd --no-check-certificate ',handles.iwg1url));
	flun3 = fopen('iwg1');
	handles.newdumpiwg = fgetl(flun3);
	fclose(flun3);
    
	if ((handles.newdumpiwg == -1))			%isempty(handles.newdumpiwg)
        display('There is currently no IWG1 data');
	else
        IWG1=textscan(handles.newdumpiwg, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', ',','EmptyValue', -Inf);   
        Time=IWG1{2};
        Latitude=IWG1{3};
        Longitude=IWG1{4};
        GPS_Alt_MSL=IWG1{5};
        GPSAltitude=IWG1{6};
        PressureAlt=IWG1{7};
        RadarAlt=IWG1{8};
        GroundSpeed=IWG1{9};
        TrueAirspeed=IWG1{10};
        IndicatedAirspeed=IWG1{11};
        Mach=IWG1{12};
        VertVelocity=IWG1{13};
        TrueHeading=IWG1{14};
        Track=IWG1{15};
        Drift=IWG1{16};
        Pitch=IWG1{17};
        Roll=IWG1{18};
        SideSlip=IWG1{19};
        AofA=IWG1{20};
        AmbientTemp=IWG1{21};
        DewPoint=IWG1{22};
        TotalTemp=IWG1{23};
        StaticPressure=IWG1{24};
        DynamicPressure=IWG1{25};
        CabinPressure=IWG1{26};
        WindSpeed=IWG1{27};
        WindDir=IWG1{28};
        VertWindSpeed=IWG1{29};
        SolarZenith=IWG1{30};
        SunElevAC=IWG1{31};
        SunAzGround=IWG1{32};
        SunAzAC=IWG1{33};
        set(handles.uitable3, 'Data',[Time; Longitude; Latitude; GPSAltitude; GPS_Alt_MSL; PressureAlt; RadarAlt; GroundSpeed;...
            TrueAirspeed; IndicatedAirspeed; Mach; VertVelocity; TrueHeading; Track; Drift; Pitch; Roll; SideSlip; AofA; AmbientTemp;...
            DewPoint; TotalTemp; StaticPressure; DynamicPressure; CabinPressure; WindSpeed; WindDir; VertWindSpeed; SolarZenith;...
            SunElevAC; SunAzGround; SunAzAC]);
        if(strcmp(handles.newTimeiwg,Time))
            set(handles.uitable3,'ForegroundColor',[0.31, 0.31, 0.31]);
        else
        set(handles.uitable3,'ForegroundColor',[0, 0, 0]);
        end;
        handles.newTimeiwg=Time;
	end;
end;         

        
    %% 
    %CPL Update
%    [handles.tempdumpcpl,handles.statuscpl]=urlreadbin('http://asp-interface-2.arc.nasa.gov/API/binary_packet_data/N806NA/CPL_BIN');
 	system(horzcat(handles.wgetpath,' wget -q -O ',handles.cplpath,' -r -nd --no-check-certificate ',handles.cplurl));
	flun2 = fopen('cpl');
	handles.tempdumpcpl = fread(flun2,1678,'uint8');
	fclose(flun2);

	   if isempty(handles.tempdumpcpl) 
            display('There is currently no CPL data');
            pause(2);
        else
            if get(handles.checkbox3, 'Value')
        %handles.newdumpcpl = uint16(handles.tempdumpcpl);
        handles.newdumpcpl = handles.tempdumpcpl';
        handles.cpltime=handles.tempdumpcpl(1:12);
        handles.cplprof(2:360,:)=handles.cplprof(1:359,:);
        for j=14:2:1678
            %handles.cpl(1,(j-12)/2)=double(handles.newdumpcpl(1,j-1))+double(bitshift(handles.newdumpcpl(1,j),7));
            handles.cpl(1,(j-12)/2)=double(handles.newdumpcpl(1,j-1))+double(handles.newdumpcpl(1,j))*255;
        end;
        handles.cplprof(1,:)=handles.cpl;
        axes(handles.axes2);
        if get(handles.radiobutton2,'Value')
            rotate3d on;
            pause(5);
        else
            rotate3d off;
            
            % set(hObject, 'String', {'Jet','HSV','Hot', 'Cool','Spring','Summer','Autumn','Winter'...
            % 'Gray','Bone','Copper','Pink','Rainbow+White'});
            clr=get(handles.popupmenu4,'Value');
            
            if (clr > 1)
            colormap(lower(handles.colors{clr}));
            else
            colormap(handles.rainbow)
            end;
            
            
            % set(hObject, 'String', {'Image', 'Surf','Mesh','Mesh Contour (meshc)', 'Mesh(meshz)','Waterfall',...
            %     'Contour Fill (contourf)','Contour 3D (contour3)', 'Surf Contour (surfc)','Ribbon' });
            switch get(handles.popupmenu3, 'Value');
                case 1
                    imagesc(handles.cplprof',[0 50]);
                case 2
                    surf(handles.cplprof);
                case 3
                    mesh(handles.cplprof);
                case 4
                    meshc(handles.cplprof);
                case 5
                    meshz(handles.cplprof);
                case 6
                    contourf(handles.cplprof);
                case 7
                    contour3(handles.cplprof);
                case 8
                    surfc(handles.cplprof);
                case 9
                    ribbon(handles.cplprof);
                case 10
                    pcolor(handles.cplprof);
            end;
            if get(handles.radiobutton1,'Value')
                axes(handles.axes2);
                if strcmp(get(colorbar,'Visible'),'off');
                    set(colorbar,'Visible','on');
                end;
            else
                axes(handles.axes2);
                colorbar('hide');
            end
            axes(handles.axes1);
            popup_sel_index = get(handles.popupmenu1, 'Value');
            if (popup_sel_index >=5)
                set(handles.edit1, 'String', GPSAltitude);
            end;
            handles.cplprof(1,:) = handles.cplprof(1,:) - mean(handles.cplprof(1,800:833));             %ADD ANOTHER RADIO BUTTON TO TOGGLE BACKGROUND SUB
            switch popup_sel_index      %'Raw', 'Background Subtracted','Raw - Altitude', 'Background Subtracted - Altitude'
                case 1
                c = plot(handles.cpl); xlim([0 833]); ylim([0 400]);
                set(c,'YData',handles.cpl);
                    pause(1);
                case 2
                    handles.cpl(1,:) = handles.cpl(1,:) - mean(handles.cpl(1,800:833)); 
                    plot(handles.cpl); xlim([0 833]); ylim([0 400]);
                    pause(1);
                case 3
                    if strcmp('NaN',num2str(str2double(get(handles.edit1,'String'))));
                        display('The Altitude must be numerical')
                        pause(1);
                    else
                        plot(handles.cpl); xlim([0 833]); ylim([0 400]);
                        set(gca,'XTickLabel',{str2double(get(handles.edit1,'String'))-((0:100:800)*30)});
                        pause(1);
                    end;
                case 4
                    if strcmp('NaN',num2str(str2double(get(handles.edit1,'String'))));
                        display('The Altitude must be numerical')
                        pause(1);
                    else
                        beep
                        pause(1);
                    end;
                case 5
                        plot(handles.cpl); xlim([0 833]); ylim([0 400]);
                        set(gca,'XTickLabel',{str2double(get(handles.edit1,'String'))-((0:100:800)*30)});
                        pause(1);
                case 6
                        handles.cpl(1,:) = handles.cpl(1,:) - mean(handles.cpl(1,800:833));         %UPDATED
                        plot(handles.cpl); xlim([0 833]); ylim([0 400]);
                        set(gca,'XTickLabel',{str2double(get(handles.edit1,'String'))-((0:100:800)*30)});
                        pause(1);
                case 7
                    pause(1);
            end;
        end;
            else
                pause(1);
            end;
        end;
    end;
else
    set(handles.pushbutton1,'String', 'Start');
end;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(handles.popupmenu1, 'Value')>=5)
    set(handles.edit1, 'Enable', 'off');
else
    set(handles.edit1, 'Enable', 'on');
end;
% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1




% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Raw', 'Background Subtracted','Raw - Manual Altitude',...
    'Background Subtracted - Manual Altitude','Raw - IWG1 Altitude', 'Background Subtracted - IWG1 Altitude', 'CAMAL Profiles' });



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    set(handles.axes1, 'Visible','on');
    set(handles.popupmenu1, 'Visible','on');
    set(handles.edit1, 'Visible','on');
    set(handles.axes2, 'Visible','on');
    set(handles.uipanel5, 'Visible','on');
    set(handles.uitable2, 'Visible','on');
else
    set(handles.axes1, 'Visible','off');
    set(handles.popupmenu1, 'Visible','off');
    set(handles.edit1, 'Visible','off');
    set(handles.axes2, 'Visible','off');
    set(handles.uipanel5, 'Visible','off');
    axes(handles.axes2);
    set(handles.uitable2, 'Visible','off');
    cla;
end

            
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    set(handles.uitable3, 'Visible','on');
else
    set(handles.uitable3, 'Visible','off');
end
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3  (CPL Checkbox).
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    set(handles.axes1, 'Visible','on');
    set(handles.popupmenu1, 'Visible','on');
    set(handles.edit1, 'Visible','on');
    set(handles.axes2, 'Visible','on');
    set(handles.uipanel5, 'Visible','on');
else
    set(handles.axes1, 'Visible','off');
    set(handles.popupmenu1, 'Visible','off');
    set(handles.edit1, 'Visible','off');
    set(handles.axes2, 'Visible','off');
    set(handles.uipanel5, 'Visible','off');
    axes(handles.axes2);
    cla;
end
% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in togglebutton1  (CAMAL or CPL Profile Toggle).
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axes(handles.axes1);
    cla;
if get(hObject,'Value')
    set(handles.togglebutton1, 'String','CPL');
%     set(handles.popupmenu1, 'Visible','off');
%     set(handles.edit1, 'Visible','off');
    set(handles.popupmenu1, 'Value', 6);            %was 7

else
    set(handles.togglebutton1, 'String','CAMAL');
     set(handles.popupmenu1, 'Visible','on');
     set(handles.edit1, 'Visible','on');
    set(handles.popupmenu1, 'Value', 1);
end
% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on selection change in popupmenu3 (Curtain Plot Type).
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Image', 'Surf','Mesh','Mesh Contour (meshc)', 'Mesh(meshz)',...
    'Contour Fill (contourf)','Contour 3D (contour3)', 'Surf Contour (surfc)','Ribbon', 'Pseudo Color (pcolor)'});


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
% (Colormap of curtain plot)
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Rainbow','Jet','HSV','Hot', 'Cool','Spring','Summer','Autumn','Winter'...
    'Gray','Bone','Copper','Pink'});


% --- Executes on button press in radiobutton1. (Curtain Plot Color Bar)
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2. (Curtain Plot rotate 3d)
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
