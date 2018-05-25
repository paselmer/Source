 handles.wgetpath = 'C:\Users\akupchoc.NDC\Documents\MobaXterm\slash\bin\wget.exe';
 handles.currentdir = cd;
 handles.camalpath = horzcat(handles.currentdir,filesep,'camal');
 handles.camalurl='https://asp-interface.arc.nasa.gov/API/binary_packet_data/N806NA/CAMAL_INMARSAT';
 handles.camalprof=zeros(360,1000);
 handles.currentdirectory = cd;
 zero4rainbow = zeros(1,21);
 one4rainbow = ones(1,21);
 taper4rainbow = (0:.05:1);
 handles.UnixStartTime = [1970 1 1 0 0 0];
 handles.GPSTime = [1980 1 6 0 0 0];
 rainbow = [zero4rainbow' zero4rainbow' taper4rainbow'; zero4rainbow' taper4rainbow' one4rainbow'; zero4rainbow' one4rainbow' 1-taper4rainbow'; taper4rainbow' one4rainbow' zero4rainbow'; one4rainbow' 1-taper4rainbow' zero4rainbow'; 1 1 1];
 handles.rainbow = rainbow;
 
 
 
            system(horzcat(handles.wgetpath,' wget -q -O ',handles.camalpath,' -r -nd --no-check-certificate ',handles.camalurl));        %not tested with CAMAL 11/23/2017
            
            flun = fopen('camal');
            handles.newdumpcamal  = fread(flun,2287,'uint8');                       %269hk + 18 + 2*bins       
            %handles.newdumpcamal  = fread(flun,287,'uint8');                       %269hk + 18       
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

            camalpacketlength = sprintf('CAMAL packet length is %i bytes long',length(handles.newdumpcamal));
            display(camalpacketlength);

               CAMAL=num2cell([APID; SequenceCount; ByteCount; datenum(UnixTime); SubSeconds; GoodCommandCount; BadCommandCount; MCSStatus; MCSState; NumMCSPacketsCollected; ControlState; ControlSubmode; ...
                 ScannerState; ScannerPos; ScannerTemps(1,1); ScannerTemps(1,2); ScannerTemps(1,3); ScannerTemps(1,4); ScannerTemps(1,5); ScannerTemps(1,6); ScannerTemps(1,7); ScannerTemps(1,8); ScannerTemps(1,9); ...
                 ScannerTemps(1,10); ScannerTemps(1,11); ScannerTemps(1,12); ScannerTemps(1,13); ScannerTemps(1,14); ScannerTemps(1,15); ScannerTemps(1,16);  IOinputs; IOoutputs; ...
                 ADvoltages(1,1); ADvoltages(1,2); ADvoltages(1,3); ADvoltages(1,4); ADvoltages(1,5); ADvoltages(1,6); ADvoltages(1,7); ADvoltages(1,8); ADvoltages(1,9); ADvoltages(1,10); ADvoltages(1,11); ...
                 ADvoltages(1,12); ADvoltages(1,13); ADvoltages(1,14); ADvoltages(1,15); ADvoltages(1,16); ADvoltages(1,17); ADvoltages(1,18); ADvoltages(1,19); ADvoltages(1,20); ADvoltages(1,21); ...
                 ADvoltages(1,22); ADvoltages(1,23); ADvoltages(1,24); ADvoltages(1,25); ADvoltages(1,26); ADvoltages(1,27); ADvoltages(1,28); ADvoltages(1,29); ADvoltages(1,30); ADvoltages(1,31); ADvoltages(1,32); ...
                 GpsWeek; GpsMilliseconds; Latitude; Longitude; Altitude; NumGPSbytesread; NASDATflags; NumIWG1pktsrec; Flags; datenum(Time);]);
%             
            switch MCSState
                case 1
                    CAMAL{9} = 'NOT_OPENED';
                case 2
                    CAMAL{9} = 'IDLE';
                case 3
                    CAMAL{9} = 'START_COLLECTING';
                case 4
                    CAMAL{9} = 'COLLECTING';
                case 5
                    CAMAL{9} = 'STOP_COLLECTING';
            end;            
            
            switch ControlState
                case 1
                    CAMAL{11} = 'INITIALIZING';
                case 2
                    CAMAL{11} = 'IDLE';
                case 3
                    CAMAL{11} = 'START_SYSTEM';
                case 4
                    CAMAL{11} = 'WAIT_FOR_LASER_ENABLE';
                case 5
                    CAMAL{11} = 'MOVE_TO_NEXT_SCANNER_POSITION';
                case 6
                    CAMAL{11} = 'COLLECT_DATA_AT_POSITION';
                case 7    
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

            handles.newTimeacamal=num2str(GpsMilliseconds);
                colormap(handles.rainbow)

            	plot(newcamalprof);
                xlabel('Bins');
                ylabel('Counts');
            end;

            
            
          