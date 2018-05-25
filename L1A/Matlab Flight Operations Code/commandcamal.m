

url = 'https://asp-interface.arc.nasa.gov/API/packet_submission/N806NA/CATS_BIN';
user = 'camal';
%packet = 'noop';
password = 'Misty6787^!'; 

url2 = 'https://asp-interface.arc.nasa.gov/API/packet_submission/N806NA/CAMAL_INMARSAT';
%[str]=urlread(url2,'post',{'user',user,'password',password,'packet',packet});
[str]=urlread(url2,'post',{'user','camal','password','Misty6787^!','packet',packet});
display(datestr(now))
display(str);

handles.wgetpath = 'C:\Users\akupchoc.NDC\Documents\MobaXterm\slash\bin\wget.exe';
system(horzcat(handles.wgetpath,' wget "',url2,'?user=camal&password=Misty6787&!" --post-data="noop" --no-check-certificate '));        %not tested with CAMAL 11/23/2017

%system('curl -k https://asp-interface.arc.nasa.gov/API/packet_submission/N806NA/CAMAL_INMARSAT -X POST -data {"user":"camal","password":"Misty6787&!","packet":"noop"} '); 


wget "https://asp-interface.arc.nasa.gov/API/packet_submission/N806NA/CAMAL_INMARSAT user=camal&password=Misty6787&! --post-data="noop" --no-check-certificate



%this worked
%wget -O test "https://asp-interface.arc.nasa.gov/API/packet_submission/N806NA/CAMAL_INMARSAT?user=camal&password=Misty6787^!" --post-data="packet=noop" --no-check-certificate

system('C:\Users\akupchoc.NDC\Documents\MobaXterm\slash\bin\wget.exe wget -O test "https://asp-interface.arc.nasa.gov/API/packet_submission/N806NA/CAMAL_INMARSAT?user=camal&password=Misty6787^!" --post-data="packet=noop" --no-check-certificate')

