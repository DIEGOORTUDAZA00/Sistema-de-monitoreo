
% Nombre del archivo: N.m
% Descripción: INTERFAZ HMI 

%% ARCHIVO Nª20 |->|
%  _______   _|_    _______
% |_______| |___|  |_______|
% Mes de Creación: MAYO
% Última Modificación: 15/09/24
% Líneas de Código: 1900



function varargout = N(varargin)


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @N_OpeningFcn, ...
                   'gui_OutputFcn',  @N_OutputFcn, ...
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


function N_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
clc;     
clear all; 
global a;  
global stop;
a= arduino('COM4','ESP32-WROOM-DevKitV1','BaudRate',9600);


function varargout = N_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function cmd_monitorear_Callback(hObject, eventdata, handles)
global stop;
stop = false;
sensor_type= get(handles.pm_sensor,'Value');
medicion_type= get(handles.pm_signal,'Value');
filtro_type= get(handles.pm_filtro,'Value');
eje_x=get(handles.cb_x,'Value');
eje_y=get(handles.cb_y,'Value');
eje_z=get(handles.cb_z,'Value');

% *************************************************************************
% * ACELEROMETRO 1
% *************************************************************************
if sensor_type==1 
   if medicion_type==1      
       global k a  
       x = 0; y=0; z=0;  
       datosx=[];  datosy=[]; datosz=[]; 
       capx=[]; capy=[]; capz=[];
       timecap=[]; 
       valor_px = 0;
       valor_ppx = 0;
       valcor_rmsx = 0;
       vel_x = 0; vel_y = 0; vel_z = 0;
       desp_x = 0; desp_y = 0; desp_z = 0;
       acel_x= 0; acel_y = 0; acel_z = 0;
       fbx=0;fby=0;fbz=0;
       timer=0;
       contador=1;
       f = 200; T = 1/f; 
       time=linspace(0, 300,500);
       waveform_x = sin(6 * pi * f .* time/T + deg2rad(45));
       waveform_y = sin(6 * pi * f .* time/T + deg2rad(90));
       waveform_z = sin(6 * pi * f .* time/T + deg2rad(180));
       hold off;
       tic;
 while toc< handles.xSamples
       if stop
           break;
       end 

       b = a.readVoltage('D36');   b=(b-1.5738)/0.33; b = b + waveform_x; 
       c = a.readVoltage('D39');   c=(c-1.5738)/0.33; c = c + waveform_y;
       d = a.readVoltage('D34');   d=(d-1.94+0.33)/0.33; d = d + waveform_z;  
       x=[x,b];  y=[y,c];      z=[z,d];        
       if filtro_type==1
       [q, t] = butter(2, 0.2, 'low');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end

       if filtro_type==2
       [q, t] = butter(2, 0.67, 'high');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end 

       if filtro_type==3
       [q, t] = butter(2,[0.2  0.9] , 'bandpass');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       hold off;
       end

       tiempo(contador) = toc;
       contador= contador+1;
              
       axes(handles.axes1);
       zoom on;
       if eje_x==1
       plot(tiempo(1:contador-1),fbx(1:contador-1),'-',Color='g'); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');  
       datosx = [datosx, b];
       acel_x = sum(abs(x));
       acel_x = sprintf('%.2f',acel_x);
       valor_px = max(abs(datosx));
       valor_px = sprintf('%.2f',valor_px);
       valor_ppx = 2*max(abs(datosx));
       valor_ppx = sprintf('%.2f',valor_ppx);
       valor_rmsx = max(abs(datosx))/sqrt(2);
       valor_rmsx = sprintf('%.2f',valor_rmsx);

       set(handles.txt_picox, 'String', valor_px);
       set(handles.txt_ppx, 'String', valor_ppx);
       set(handles.txt_rmsx,'String',valor_rmsx);
       set(handles.txt_medx, 'String', acel_x);
       fxcap=fbx;
       if toc>=0 && toc<=1
       capx=[capx,fxcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picox, 'String', '');
               set(handles.txt_ppx, 'String', '');
               set(handles.txt_rmsx, 'String', '');
               set(handles.txt_medx, 'String', '');
       end

       if eje_y==1
       plot(tiempo(1:contador-1),fby(1:contador-1),'-','Color',[1, 0.5, 0],'LineWidth',0.1); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosy = [datosy, c];
       acel_y=sum(abs(y));
       acel_y = sprintf('%.2f',acel_y);
       valor_py = max(abs(datosy));
       valor_py = sprintf('%.2f',valor_py);
       valor_ppy = 2*max(abs(datosy));
       valor_ppy = sprintf('%.2f',valor_ppy);
       valor_rmsy = max(abs(datosy))/sqrt(2);
       valor_rmsy = sprintf('%.2f',valor_rmsy);
       set(handles.txt_picoy, 'String', valor_py);
       set(handles.txt_ppy, 'String', valor_ppy);
       set(handles.txt_rmsy,'String',valor_rmsy);
       set(handles.txt_medy, 'String', acel_y);
       fycap=fby;      
       if toc>=0 && toc<=1  
       capy=[capy,fycap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picoy, 'String', '');
               set(handles.txt_ppy, 'String', '');
               set(handles.txt_rmsy, 'String', '');
               set(handles.txt_medy, 'String', '');
       end
     
       if eje_z==1
       plot(tiempo(1:contador-1),fbz(1:contador-1),Color='b');hold on;       
       xlabel(handles.axes1,'FRECUENCIA' );
       ylabel(handles.axes1,'AMPLITUD' );
       datosz = [datosz, d];       
       acel_z = sum(abs(z));
       acel_z = sprintf('%.2f',acel_z);
       valor_pz = max(abs(datosz));
       valor_pz = sprintf('%.2f',valor_pz);
       valor_ppz = 2*max(abs(datosz));
       valor_ppz = sprintf('%.2f',valor_ppz);
       valor_rmsz = max(abs(datosz))/sqrt(2);
       valor_rmsz = sprintf('%.2f',valor_rmsz);
       set(handles.txt_picoz, 'String', valor_pz);
       set(handles.txt_ppz, 'String', valor_ppz);
       set(handles.txt_rmsz,'String',valor_rmsz);
       set(handles.txt_medz, 'String', acel_z);      
       fzcap=fbz;
       if toc>=0 && toc<=1
       capz=[capz,fzcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picoz, 'String', '');
               set(handles.txt_ppz, 'String', '');
               set(handles.txt_rmsz, 'String', '');
               set(handles.txt_medz, 'String', '');
       end
       axis([1 handles.xSamples -3 3]);
       drawnow;
     save('senalx',"capx");  
     save('senaly',"capy");
     save('senalz',"capz");  
   t=8;
   Lx=numel(capx);
   Ly=numel(capy);
   Lz=numel(capz);
   Fsx=Lx/t;
   Fsy=Ly/t;
   Fsz=Lz/t;
   fx=Fsx/Lx*(0:Lx-1);
   fy=Fsy/Ly*(0:Ly-1);
   fz=Fsz/Lz*(0:Lz-1);
   ftx = fft(capx); 
   fty = fft(capy);
   ftz = fft(capz);
        axes(handles.axes2);

        if eje_x==1   
        plot(tiempo(1:contador-1),x(1:contador-1),Color='g'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else     
        end 
        if eje_y==1   
        plot(tiempo(1:contador-1),y(1:contador-1),'Color',[1, 0.5, 0]); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else  
        end
        if eje_z==1
        plot(tiempo(1:contador-1),z(1:contador-1),'b'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else
        end
       axis([1 handles.xSamples -3 3]);    
       pause(1);
 end

if eje_x == 1 
figure;
plot(fx,abs(ftx),Color='g');
xlim([1 88]);
saveas(gcf, 'espectro_x.jpg');
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje X');  
end

if eje_y == 1
figure;
plot(fy, abs(fty), 'Color', 'r'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Y');
saveas(gcf, 'espectro_y.jpg');
end

if eje_z == 1
figure;
plot(fz, abs(ftz), 'Color', 'b');  
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Z');
saveas(gcf, 'espectro_z.jpg');
end
end
   
   if medicion_type==2
       global k a 
       x = 0; y=0; z=0;    
       datosx=[];  datosy=[]; datosz=[]; 
       capx=[]; capy=[]; capz=[];
       timecap=[]; 
       valor_px = 0;
       valor_ppx = 0;
       valor_rmsx = 0;
       vel_x = 0; vel_y = 0; vel_z = 0;
       desp_x = 0; desp_y = 0; desp_z = 0;
       acel_x= 0; acel_y = 0; acel_z = 0;
       fbx=0;fby=0;fbz=0;
       timer=0;
       contador=1;
       f = 36.8; T = 1/f; 
       time=linspace(0, 300,500);
       waveform_x = sin(2 * pi * f .* time/T + deg2rad(45));
       waveform_y = sin(2 * pi * f .* time/T + deg2rad(90));
       waveform_z = sin(2 * pi * f .* time/T + deg2rad(180));
       hold off;
       tic;
   while toc < handles.xSamples
       if stop
           break;
       end
       b = a.readVoltage('D36');   b=(b-1.5738)/0.33;  b = b + waveform_x;   
       c = a.readVoltage('D39');   c=(c-1.5738)/0.33;  c = c + waveform_y;
       d = a.readVoltage('D34');   d=(d-1.94+0.33)/0.33;  d = d + waveform_z;   
       x=[x,b];  y=[y,c];      z=[z,d];     

       if filtro_type==1
       [q, t] = butter(2, 0.2, 'low');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end

       if filtro_type==2
       [q, t] = butter(2, 0.67, 'high');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end 

       if filtro_type==3
       [q, t] = butter(2,[0.2  0.9] , 'bandpass');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       hold off;
       end

       tiempo(contador) = toc;
       contador= contador+1; 

       axes(handles.axes1);
       zoom on;
       if eje_x==1  
       plot(tiempo(1:contador-1),fbx(1:contador-1),'-',Color='g'); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosx = [datosx, b];
       acel_x = sum(abs(x));
       acel_x = sprintf('%.2f',acel_x);
       valor_px = max(abs(datosx));
       valor_px = (1000*9.8*valor_px)/(2*pi*1000); 
       valor_ppx = 2*max(abs(valor_px));
       valor_rmsx = (valor_px)/(sqrt(2));
       valor_ppx = sprintf('%.2f',valor_ppx);
       valor_rmsx = sprintf('%.2f',valor_rmsx);
       valor_px = sprintf('%.2f',valor_px);
       set(handles.txt_picox, 'String', valor_px);
       set(handles.txt_ppx, 'String', valor_ppx);
       set(handles.txt_rmsx,'String',valor_rmsx);
       set(handles.txt_medx, 'String', acel_x);
       fxcap=fbx;
       if toc>=0 && toc<=1
       capx=[capx,fxcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end

       else
               set(handles.txt_picox, 'String', '');
               set(handles.txt_ppx, 'String', '');
               set(handles.txt_rmsx, 'String', '');
               set(handles.txt_medx, 'String', '');
       end
        
       if eje_y==1
       plot(tiempo(1:contador-1),fby(1:contador-1),'-','Color',[1, 0.5, 0],'LineWidth',0.1); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosy = [datosy, c];
       acel_y=sum(abs(y));
       acel_y = sprintf('%.2f',acel_y);
       valor_py = max(abs(datosy));
       valor_py = (1000*9.8*valor_py)/(2*pi*1000);
       valor_ppy = 2*max(abs(datosy));
       valor_rmsy = (valor_py)/(sqrt(2));
       valor_rmsy = sprintf('%.2f',valor_rmsy);
       valor_py = sprintf('%.2f',valor_py);
       valor_ppy = sprintf('%.2f',valor_ppy);
       set(handles.txt_picoy, 'String', valor_py);
       set(handles.txt_ppy, 'String', valor_ppy);
       set(handles.txt_rmsy,'String',valor_rmsy);
       set(handles.txt_medy, 'String', acel_y);
       fycap=fby;
       if toc>=0 && toc<=1
           capy=[capy,fycap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picoy, 'String', '');
               set(handles.txt_ppy, 'String', '');
               set(handles.txt_rmsy, 'String', '');
               set(handles.txt_medy, 'String', '');
       end
     
       if eje_z==1
       plot(tiempo(1:contador-1),fbz(1:contador-1),Color='b');hold on;
       xlabel(handles.axes1,'FRECUENCIA' );
       ylabel(handles.axes1,'AMPLITUD' );
       datosz = [datosz, d];
       acel_z = sum(abs(z));
       acel_z = sprintf('%.2f',acel_z);
       valor_pz = max(abs(datosz));
       valor_pz = (1000*9.8*valor_pz)/(2*pi*1000);
       valor_ppz = 2*max(abs(datosz));
       valor_rmsz = (valor_pz)/(sqrt(2));
       valor_pz = sprintf('%.2f',valor_pz);
       valor_ppz = sprintf('%.2f',valor_ppz);
       valor_rmsz = sprintf('%.2f',valor_rmsz);
       set(handles.txt_picoz, 'String', valor_pz);
       set(handles.txt_ppz, 'String', valor_ppz);
       set(handles.txt_rmsz,'String',valor_rmsz);
       set(handles.txt_medz, 'String', acel_z);
       fzcap=fbz;
       if toc>=0 && toc<=1
           capz=[capz,fzcap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picoz, 'String', '');
               set(handles.txt_ppz, 'String', '');
               set(handles.txt_rmsz, 'String', '');
               set(handles.txt_medz, 'String', '');
       end
       axis([1 handles.xSamples -3 3]);
       drawnow;
       save('senalx',"capx");  
       save('senaly',"capy");
       save('senalz',"capz");
   t=8;
   Lx=numel(capx);
   Ly=numel(capy);
   Lz=numel(capz);
   Fsx=Lx/t;
   Fsy=Ly/t;
   Fsz=Lz/t;
   fx=Fsx/Lx*(0:Lx-1);
   fy=Fsy/Ly*(0:Ly-1);
   fz=Fsz/Lz*(0:Lz-1);
   ftx = fft(capx); 
   fty = fft(capy);
   ftz = fft(capz);

       axes(handles.axes2);
        if eje_x==1   
        plot(tiempo(1:contador-1),x(1:contador-1),Color='g'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else 
          
           
        end 
        if eje_y==1
        plot(tiempo(1:contador-1),y(1:contador-1),'Color',[1, 0.5, 0]); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else
            
        end
        if eje_z==1
        plot(tiempo(1:contador-1),z(1:contador-1),'b'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else
         
        end
       axis([1 handles.xSamples -3 3]);               
       pause(1); 
   end   
   if eje_x == 1 
figure;
plot(fx,abs(ftx),Color='g');
xlim([1 88]);
saveas(gcf, 'espectro_x.jpg');
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje X');  
end

if eje_y == 1
figure;
plot(fy, abs(fty), 'Color', 'r'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Y');
saveas(gcf, 'espectro_y.jpg');
end

if eje_z == 1
figure;
plot(fz, abs(ftz), 'Color', 'b');  
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Z');
saveas(gcf, 'espectro_z.jpg');
end
   end
   if medicion_type==3
       global k a  
       x = 0; y=0; z=0;
       datosx=[];  datosy=[]; datosz=[]; 
       capx=[]; capy=[]; capz=[];
       timecap=[];
       valor_px = 0;
       valor_ppx = 0;
       valor_rmsx = 0;
       vel_x = 0; vel_y = 0; vel_z = 0;
       desp_x = 0; desp_y = 0; desp_z = 0;
       acel_x= 0; acel_y = 0; acel_z = 0;
       fbx=0;fby=0;fbz=0;
       contador=1;
       f = 200; T = 1/f; 
       time=linspace(0, 300,500);
       waveform_x = sin(6 * pi * f .* time/T + deg2rad(45));
       waveform_y = sin(6 * pi * f .* time/T + deg2rad(90));
       waveform_z = sin(6 * pi * f .* time/T + deg2rad(180));
       hold off;
       tic;
       while toc < handles.xSamples
 
       if stop
           break;
       end
       b = a.readVoltage('D36');   b=(b-1.5738)/0.33;   b = b + waveform_x;  
       c = a.readVoltage('D39');   c=(c-1.5738)/0.33;   c = c + waveform_y;
       d = a.readVoltage('D34');   d=(d-1.94+0.33)/0.33;  d = d + waveform_z;  
       x=[x,b];  y=[y,c];      z=[z,d];      

       if filtro_type==1
       [q, t] = butter(2, 0.2, 'low');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end

       if filtro_type==2
       [q, t] = butter(2, 0.67, 'high');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end 

       if filtro_type==3
       [q, t] = butter(2,[0.2  0.9] , 'bandpass');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       hold off;
       end
       
       tiempo(contador) = toc;
       contador= contador+1;

       axes(handles.axes1);
       zoom on;

       if eje_x==1  
       plot(tiempo(1:contador-1),fbx(1:contador-1),'-',Color='g'); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosx = [datosx, b];
       acel_x = sum(abs(x));
       acel_x = sprintf('%.2f',acel_x);
       valor_px = max(abs(datosx));
       valor_ppx = valor_px/(pi*10);
       valor_rmsx = (valor_ppx)/(2*sqrt(2));
       valor_ppx = sprintf('%.2f',valor_ppx);
       valor_rmsx = sprintf('%.2f',valor_rmsx);
       valor_px = sprintf('%.2f',valor_px);
       set(handles.txt_picox, 'String', valor_px);
       set(handles.txt_ppx, 'String', valor_ppx);
       set(handles.txt_rmsx,'String',valor_rmsx);
       set(handles.txt_medx, 'String', acel_x);
       fxcap=fbx;
       if toc>=0 && toc<=1
              capx=[capx,fxcap];
              set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picox, 'String', '');
               set(handles.txt_ppx, 'String', '');
               set(handles.txt_rmsx, 'String', '');
               set(handles.txt_medx, 'String', '');
       end
        
       if eje_y==1
       plot(tiempo(1:contador-1),fby(1:contador-1),'-','Color',[1, 0.5, 0],'LineWidth',0.1); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosy = [datosy, c];
       acel_y=sum(abs(y));
       acel_y = sprintf('%.2f',acel_y);
       valor_py = max(abs(datosy));
       valor_ppy = (valor_py)/(pi*10);
       valor_rmsy = (valor_ppy)/(2*sqrt(2));
       valor_rmsy = sprintf('%.2f',valor_rmsy);
       valor_py = sprintf('%.2f',valor_py);
       valor_ppy = sprintf('%.2f',valor_ppy);
       set(handles.txt_picoy, 'String', valor_py);
       set(handles.txt_ppy, 'String', valor_ppy);
       set(handles.txt_rmsy,'String',valor_rmsy);
       set(handles.txt_medy, 'String', acel_y);
       fycap=fby;
       if toc>=0 && toc<=1
           capy=[capy,fycap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoy, 'String', '');
               set(handles.txt_ppy, 'String', '');
               set(handles.txt_rmsy, 'String', '');
               set(handles.txt_medy, 'String', '');
       end
     
       if eje_z==1
       plot(tiempo(1:contador-1),fbz(1:contador-1),Color='b');hold on;
       xlabel(handles.axes1,'FRECUENCIA' );
       ylabel(handles.axes1,'AMPLITUD' );
       datosz = [datosz, d];
       acel_z = sum(abs(z));
       acel_z = sprintf('%.2f',acel_z);
       valor_pz = max(abs(datosz));
       valor_ppz = (valor_pz)/(pi*10);
       valor_rmsz = (valor_ppz)/(2*sqrt(2));
       valor_pz = sprintf('%.2f',valor_pz);
       valor_ppz = sprintf('%.2f',valor_ppz);
       valor_rmsz = sprintf('%.2f',valor_rmsz);
       set(handles.txt_picoz, 'String', valor_pz);
       set(handles.txt_ppz, 'String', valor_ppz);
       set(handles.txt_rmsz,'String',valor_rmsz);
       set(handles.txt_medz, 'String', acel_z);
       fzcap=fbz;
       if toc>=0 && toc<=1
       capz=[capz,fzcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoz, 'String', '');
               set(handles.txt_ppz, 'String', '');
               set(handles.txt_rmsz, 'String', '');
               set(handles.txt_medz, 'String', '');
       end
       axis([1 handles.xSamples -2 2]);
       drawnow;
       save('senalx',"capx");  
       save('senaly',"capy");
       save('senalz',"capz");
       t=8;
       Lx=numel(capx);
       Ly=numel(capy);
       Lz=numel(capz);
       Fsx=Lx/t;
       Fsy=Ly/t;
       Fsz=Lz/t;
       fx=Fsx/Lx*(0:Lx-1);
       fy=Fsy/Ly*(0:Ly-1);
       fz=Fsz/Lz*(0:Lz-1);
       ftx = fft(capx); 
       fty = fft(capy);
       ftz = fft(capz);

       axes(handles.axes2);
       
        if eje_x==1 
        plot(tiempo(1:contador-1),x(1:contador-1),Color='g'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else            
        end 
        if eje_y==1
        plot(tiempo(1:contador-1),y(1:contador-1),'Color',[1, 0.5, 0]); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');             
        else            
        end
        if eje_z==1
        plot(tiempo(1:contador-1),z(1:contador-1),'b'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');      
        else         
        end
       axis([1 handles.xSamples -3 3]);  
       pause(1); 
         end 
if eje_x == 1 
figure;
plot(fx,abs(ftx),Color='g');
xlim([1 88]);
saveas(gcf, 'espectro_x.jpg');
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje X');  
end

if eje_y == 1
figure;
plot(fy, abs(fty), 'Color', 'r'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Y');
saveas(gcf, 'espectro_y.jpg');
end

if eje_z == 1
figure;
plot(fz, abs(ftz), 'Color', 'b');
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Z');
saveas(gcf, 'espectro_z.jpg');
end
   end
end
% *************************************************************************
% * ACELEROMETRO 2
% ************************************************************************* 
if sensor_type==2 
if medicion_type==1 
   global k a  
   x = 0; y=0; z=0; 
   datosx=[];  datosy=[]; datosz=[]; 
   capx=[]; capy=[]; capz=[];
   timecap=[]; 
       valor_px = 0;
       valor_ppx = 0;
       valor_rmsx = 0;
       delta_t = 1/handles.xSamples;
       vel_x = 0; vel_y = 0; vel_z = 0;
       desp_x = 0; desp_y = 0; desp_z = 0;
       acel_x= 0; acel_y = 0; acel_z = 0;
       fbx=0;fby=0;fbz=0;
       timer=0;
       contador=1;
       f = 36.8; T = 1/f; 
       time=linspace(0, 300,500);
       waveform_x = sin(2 * pi * f .* time/T + deg2rad(45));
       waveform_y = sin(2 * pi * f .* time/T + deg2rad(90));
       waveform_z = sin(2 * pi * f .* time/T + deg2rad(180));
       hold off;
       tic;
       while toc < handles.xSamples 

       if stop
           break;
       end

       b = a.readVoltage('D35');   b=(b-1.5738)/0.33;    b = b + waveform_x;
       c = a.readVoltage('D32');   c=(c-1.5738)/0.33;    c = c + waveform_y;
       d = a.readVoltage('D33');   d=(d-1.94+0.33)/0.33;  d = d + waveform_z;  
       x=[x,b];  y=[y,c];      z=[z,d];                 
       if filtro_type==1
       [q, t] = butter(2, 0.2, 'low');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end

       if filtro_type==2
       [q, t] = butter(2, 0.67, 'high');
       fbx=filter(q,t,x);%%
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end 

       if filtro_type==3
       [q, t] = butter(2,[0.2  0.9] , 'bandpass');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       hold off;
       end

       tiempo(contador) = toc;
       contador= contador+1;
              
       axes(handles.axes1);
       zoom on;

       if eje_x==1 
       plot(tiempo(1:contador-1),fbx(1:contador-1),'-',Color='g'); hold on; 
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');  
       datosx = [datosx, b];
       acel_x = sum(abs(x));
       acel_x = sprintf('%.2f',acel_x);
       valor_px = max(abs(datosx));
       valor_px = sprintf('%.2f',valor_px);
       valor_ppx = 2*max(abs(datosx));
       valor_ppx = sprintf('%.2f',valor_ppx);
       valor_rmsx = max(abs(datosx))/sqrt(2);
       valor_rmsx = sprintf('%.2f',valor_rmsx);
       set(handles.txt_picox, 'String', valor_px);
       set(handles.txt_ppx, 'String', valor_ppx);
       set(handles.txt_rmsx,'String',valor_rmsx);
       set(handles.txt_medx, 'String', acel_x);
       fxcap=fbx;
       if toc>=0 && toc<=1
           capx=[capx,fxcap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picox, 'String', '');
               set(handles.txt_ppx, 'String', '');
               set(handles.txt_rmsx, 'String', '');
               set(handles.txt_medx, 'String', '');
       end
        
       if eje_y==1
       plot(tiempo(1:contador-1),fby(1:contador-1),'-','Color',[1, 0.5, 0],'LineWidth',0.1); hold on;
       xlabel(handles.axes1,'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosy = [datosy, c];
       acel_y=sum(abs(y));
       acel_y = sprintf('%.2f',acel_y);
       valor_py = max(abs(datosy));
       valor_py = sprintf('%.2f',valor_py);
       valor_ppy = 2*max(abs(datosy));
       valor_ppy = sprintf('%.2f',valor_ppy);
       valor_rmsy = max(abs(datosy))/sqrt(2);
       valor_rmsy = sprintf('%.2f',valor_rmsy);
       set(handles.txt_picoy, 'String', valor_py);
       set(handles.txt_ppy, 'String', valor_ppy);
       set(handles.txt_rmsy,'String',valor_rmsy);
       set(handles.txt_medy, 'String', acel_y);
       fycap=fby;
       if toc>=0 && toc<=1 
           capy=[capy,fycap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoy, 'String', '');
               set(handles.txt_ppy, 'String', '');
               set(handles.txt_rmsy, 'String', '');
               set(handles.txt_medy, 'String', '');
       end
     
       if eje_z==1
       plot(tiempo(1:contador-1),fbz(1:contador-1),Color='b');hold on;
       xlabel(handles.axes1, 'FRECUENCIA' );
       ylabel(handles.axes1, 'AMPLITUD');
       datosz = [datosz, d];
       acel_z = sum(abs(z));
       acel_z = sprintf('%.2f',acel_z);
       valor_pz = max(abs(datosz));
       valor_pz = sprintf('%.2f',valor_pz);
       valor_ppz = 2*max(abs(datosz));
       valor_ppz = sprintf('%.2f',valor_ppz);
       valor_rmsz = max(abs(datosz))/sqrt(2);
       valor_rmsz = sprintf('%.2f',valor_rmsz);
       set(handles.txt_picoz, 'String', valor_pz);
       set(handles.txt_ppz, 'String', valor_ppz);
       set(handles.txt_rmsz,'String',valor_rmsz);
       set(handles.txt_medz, 'String', acel_z);
       fzcap=fbz;
       if toc>=0 && toc<=1
       capz=[capz,fzcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoz, 'String', '');
               set(handles.txt_ppz, 'String', '');
               set(handles.txt_rmsz, 'String', '');
               set(handles.txt_medz, 'String', '');
       end
       axis([1 handles.xSamples -3 3]);
       drawnow;
     save('senalx',"capx");  
     save('senaly',"capy");
     save('senalz',"capz");  
   t=8;
   Lx=numel(capx);
   Ly=numel(capy);
   Lz=numel(capz);
   Fsx=Lx/t;
   Fsy=Ly/t;
   Fsz=Lz/t;
   fx=Fsx/Lx*(0:Lx-1);
   fy=Fsy/Ly*(0:Ly-1);
   fz=Fsz/Lz*(0:Lz-1);
   ftx = fft(capx); 
   fty = fft(capy);
   ftz = fft(capz);

       axes(handles.axes2);
        if eje_x==1 
        plot(tiempo(1:contador-1),x(1:contador-1),Color='g'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else 
          
           
        end 
        if eje_y==1
        plot(tiempo(1:contador-1),y(1:contador-1),'Color',[1, 0.5, 0]); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
       
        
      
        else
            
        end
        if eje_z==1
        plot(tiempo(1:contador-1),z(1:contador-1),'b'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
      
        else
         
        end
       axis([1 handles.xSamples -3 3]); 
       pause(1);      
       end
if eje_x == 1 
figure;
plot(fx,abs(ftx),Color='g');
xlim([1 88]);
saveas(gcf, 'espectro_x.jpg');
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje X');  
end

if eje_y == 1
figure;
plot(fy, abs(fty), 'Color', 'r'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Y');
saveas(gcf, 'espectro_y.jpg');
end

if eje_z == 1
figure;
plot(fz, abs(ftz), 'Color', 'b');  
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Z');
saveas(gcf, 'espectro_z.jpg');
end       
end
   if medicion_type==2
global k a  
x = 0; y=0; z=0;    
datosx=[];  datosy=[]; datosz=[]; 
capx=[]; capy=[]; capz=[];
timecap=[];
       valor_px = 0;
       valor_ppx = 0;
       valor_rmsx = 0;
       delta_t = 1/handles.xSamples;
       vel_x = 0; vel_y = 0; vel_z = 0;
       desp_x = 0; desp_y = 0; desp_z = 0;
       acel_x= 0; acel_y = 0; acel_z = 0;
       fbx=0;fby=0;fbz=0;
       timer=0;
       contador=1;
       f = 36.8; T = 1/f; 
       time=linspace(0, 300,500);
       waveform_x = sin(2 * pi * f .* time/T + deg2rad(45));
       waveform_y = sin(2 * pi * f .* time/T + deg2rad(90));
       waveform_z = sin(2 * pi * f .* time/T + deg2rad(180));
   hold off;
   tic;
   while toc < handles.xSamples

       if stop
           break;
       end
       b = a.readVoltage('D35');   b=(b-1.5738)/0.33;    b = b + waveform_x;
       c = a.readVoltage('D32');   c=(c-1.5738)/0.33;    c = c + waveform_y;
       d = a.readVoltage('D33');   d=(d-1.94+0.33)/0.33;  d = d + waveform_z;  
       x=[x,b];  y=[y,c];      z=[z,d];     
        
       if filtro_type==1
       [q, t] = butter(2, 0.2, 'low');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end

       if filtro_type==2
       [q, t] = butter(2, 0.67, 'high');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end 

       if filtro_type==3
       [q, t] = butter(2,[0.2  0.9] , 'bandpass');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       hold off;
       end

       tiempo(contador) = toc;
       contador= contador+1;

       axes(handles.axes1);
       zoom on;

       if eje_x==1   
       plot(tiempo(1:contador-1),fbx(1:contador-1),'-',Color='g'); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosx = [datosx, b];
       acel_x = sum(abs(x));
       acel_x = sprintf('%.2f',acel_x);
       valor_px = max(abs(datosx));
       valor_px = (1000*9.8*valor_px)/(2*pi*1000);
       valor_ppx = 2*max(abs(valor_px));
       valor_rmsx = (valor_px)/(sqrt(2));
       valor_ppx = sprintf('%.2f',valor_ppx);
       valor_rmsx = sprintf('%.2f',valor_rmsx);
       valor_px = sprintf('%.2f',valor_px);
       set(handles.txt_picox, 'String', valor_px);
       set(handles.txt_ppx, 'String', valor_ppx);
       set(handles.txt_rmsx,'String',valor_rmsx);
       set(handles.txt_medx, 'String', acel_x);
       fxcap=fbx;
       if toc>=0 && toc<=1
       capx=[capx,fxcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end

       else
               set(handles.txt_picox, 'String', '');
               set(handles.txt_ppx, 'String', '');
               set(handles.txt_rmsx, 'String', '');
               set(handles.txt_medx, 'String', '');
       end
        
       if eje_y==1
       plot(tiempo(1:contador-1),fby(1:contador-1),'-','Color',[1, 0.5, 0],'LineWidth',0.1); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosy = [datosy, c];
       acel_y=sum(abs(y));
       acel_y = sprintf('%.2f',acel_y);
       valor_py = max(abs(datosy));
       valor_py = (1000*9.8*valor_py)/(2*pi*1000);
       valor_ppy = 2*max(abs(datosy));
       valor_rmsy = (valor_py)/(sqrt(2));
       valor_rmsy = sprintf('%.2f',valor_rmsy);
       valor_py = sprintf('%.2f',valor_py);
       valor_ppy = sprintf('%.2f',valor_ppy);
       set(handles.txt_picoy, 'String', valor_py);
       set(handles.txt_ppy, 'String', valor_ppy);
       set(handles.txt_rmsy,'String',valor_rmsy);
       set(handles.txt_medy, 'String', acel_y);
       fycap=fby;
       if toc>=0 && toc<=1
       capy=[capy,fycap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picoy, 'String', '');
               set(handles.txt_ppy, 'String', '');
               set(handles.txt_rmsy, 'String', '');
               set(handles.txt_medy, 'String', '');
       end
     
       if eje_z==1
       plot(tiempo(1:contador-1),fbz(1:contador-1),Color='b');hold on;
       xlabel(handles.axes1,'FRECUENCIA' );
       ylabel(handles.axes1,'AMPLITUD' );
       datosz = [datosz, d];
       acel_z = sum(abs(z));
       acel_z = sprintf('%.2f',acel_z);
       valor_pz = max(abs(datosz));
       valor_pz = (1000*9.8*valor_pz)/(2*pi*1000);
       valor_ppz = 2*max(abs(datosz));
       valor_rmsz = (valor_pz)/(sqrt(2));
       valor_pz = sprintf('%.2f',valor_pz);
       valor_ppz = sprintf('%.2f',valor_ppz);
       valor_rmsz = sprintf('%.2f',valor_rmsz);
       set(handles.txt_picoz, 'String', valor_pz);
       set(handles.txt_ppz, 'String', valor_ppz);
       set(handles.txt_rmsz,'String',valor_rmsz);
       set(handles.txt_medz, 'String', acel_z);
       fzcap=fbz;
       if toc>=0 && toc<=1
           capz=[capz,fzcap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoz, 'String', '');
               set(handles.txt_ppz, 'String', '');
               set(handles.txt_rmsz, 'String', '');
               set(handles.txt_medz, 'String', '');
       end
       axis([1 handles.xSamples -3 3]);
       drawnow;
       save('senalx',"capx");  
       save('senaly',"capy");
       save('senalz',"capz");
   t=8;
   Lx=numel(capx);
   Ly=numel(capy);
   Lz=numel(capz);
   Fsx=Lx/t;
   Fsy=Ly/t;
   Fsz=Lz/t;
   fx=Fsx/Lx*(0:Lx-1);
   fy=Fsy/Ly*(0:Ly-1);
   fz=Fsz/Lz*(0:Lz-1);
   ftx = fft(capx); 
   fty = fft(capy);
   ftz = fft(capz);

       axes(handles.axes2);
       
        if eje_x==1  
        plot(tiempo(1:contador-1),x(1:contador-1),Color='g'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else 
          
           
        end 
        if eje_y==1
        plot(tiempo(1:contador-1),y(1:contador-1),'Color',[1, 0.5, 0]); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
       
        
      
        else
            
        end
        if eje_z==1
        plot(tiempo(1:contador-1),z(1:contador-1),'b'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
      
        else
         
        end
       axis([1 handles.xSamples -3 3]);   
       pause(1); 
   end    
if eje_x == 1 
figure;
plot(fx,abs(ftx),Color='g');
xlim([1 88]);
saveas(gcf, 'espectro_x.jpg');
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje X');  
end

if eje_y == 1
figure;
plot(fy, abs(fty), 'Color', 'r'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Y');
saveas(gcf, 'espectro_y.jpg');
end

if eje_z == 1
figure;
plot(fz, abs(ftz), 'Color', 'b'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Z');
saveas(gcf, 'espectro_z.jpg');
end
end
   
   if medicion_type==3
global k a 
   x = 0; y=0; z=0;    
   datosx=[];  datosy=[]; datosz=[];
   capx=[]; capy=[]; capz=[];
   timecap=[];
       valor_px = 0;
       valor_ppx = 0;
       valor_rmsx = 0;
       vel_x = 0; vel_y = 0; vel_z = 0;
       desp_x = 0; desp_y = 0; desp_z = 0;
       acel_x= 0; acel_y = 0; acel_z = 0;
       fbx=0;fby=0;fbz=0;
       timer=0;
       contador=1;
       f = 36.8; T = 1/f; 
       time=linspace(0, 300,500);
       waveform_x = sin(2 * pi * f .* time/T + deg2rad(45));
       waveform_y = sin(2 * pi * f .* time/T + deg2rad(90));
       waveform_z = sin(2 * pi * f .* time/T + deg2rad(180));
       hold off;
       tic;
       while toc < handles.xSamples 
   
       if stop
           break;
       end
       %% HASTA AQUI VIENE LA PARTE DE OBTENCION DE DATOS DEL ACELEROMETRO 1 Y 2
       b = a.readVoltage('D35');   b=(b-1.5738)/0.33;   b = b + waveform_x;  
       c = a.readVoltage('D32');   c=(c-1.5738)/0.33;   c = c + waveform_y;
       d = a.readVoltage('D33');   d=(d-1.94+0.33)/0.33;  d = d + waveform_z;  
       x=[x,b];  y=[y,c];      z=[z,d];                
      
       if filtro_type==1
       [q, t] = butter(2, 0.2, 'low');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end

       if filtro_type==2
       [q, t] = butter(2, 0.67, 'high');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       end 

       if filtro_type==3
       [q, t] = butter(2,[0.2  0.67] , 'bandpass');
       fbx=filter(q,t,x);
       fby=filter(q,t,y);
       fbz=filter(q,t,z);
       hold off;
       end

       tiempo(contador) = toc;
       contador= contador+1;

       axes(handles.axes1);
       zoom on;

       if eje_x==1  
       plot(tiempo(1:contador-1),fbx(1:contador-1),'-',Color='g'); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosx = [datosx, b];
       acel_x = sum(abs(x));
       acel_x = sprintf('%.2f',acel_x);
       valor_px = max(abs(datosx));
       valor_ppx = valor_px/(pi*10);
       valor_rmsx = (valor_ppx)/(2*sqrt(2));
       valor_ppx = sprintf('%.2f',valor_ppx);
       valor_rmsx = sprintf('%.2f',valor_rmsx);
       valor_px = sprintf('%.2f',valor_px);
       set(handles.txt_picox, 'String', valor_px);
       set(handles.txt_ppx, 'String', valor_ppx);
       set(handles.txt_rmsx,'String',valor_rmsx);
       set(handles.txt_medx, 'String', acel_x);
       fxcap=fbx;
       if toc>=0 && toc<=1
       capx=[capx,fxcap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
       else
               set(handles.txt_picox, 'String', '');
               set(handles.txt_ppx, 'String', '');
               set(handles.txt_rmsx, 'String', '');
               set(handles.txt_medx, 'String', '');
       end
        
       if eje_y==1
       plot(tiempo(1:contador-1),fby(1:contador-1),'-','Color',[1, 0.5, 0],'LineWidth',0.1); hold on;
       xlabel(handles.axes1, 'FRECUENCIA');
       ylabel(handles.axes1, 'AMPLITUD');
       datosy = [datosy, c];
       acel_y=sum(abs(y));
       acel_y = sprintf('%.2f',acel_y);
       valor_py = max(abs(datosy));
       valor_ppy = (valor_py)/(pi*10);
       valor_rmsy = (valor_ppy)/(2*sqrt(2));
       valor_rmsy = sprintf('%.2f',valor_rmsy);
       valor_py = sprintf('%.2f',valor_py);
       valor_ppy = sprintf('%.2f',valor_ppy);
       set(handles.txt_picoy, 'String', valor_py);
       set(handles.txt_ppy, 'String', valor_ppy);
       set(handles.txt_rmsy,'String',valor_rmsy);
       set(handles.txt_medy, 'String', acel_y);
       fycap=fby;
       time=contador;
       if toc>=0 && toc<=1
       capy=[capy,fycap];
       set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoy, 'String', '');
               set(handles.txt_ppy, 'String', '');
               set(handles.txt_rmsy, 'String', '');
               set(handles.txt_medy, 'String', '');
       end
     
       if eje_z==1
       plot(tiempo(1:contador-1),fbz(1:contador-1),Color='b');hold on;    
       xlabel(handles.axes1,'FRECUENCIA' );
       ylabel(handles.axes1,'AMPLITUD' );
       datosz = [datosz, d];
       acel_z = sum(abs(z));
       acel_z = sprintf('%.2f',acel_z);
       valor_pz = max(abs(datosz));
       valor_ppz = (valor_pz)/(pi*10);
       valor_rmsz = (valor_ppz)/(2*sqrt(2));
       valor_pz = sprintf('%.2f',valor_pz);
       valor_ppz = sprintf('%.2f',valor_ppz);
       valor_rmsz = sprintf('%.2f',valor_rmsz);
       set(handles.txt_picoz, 'String', valor_pz);
       set(handles.txt_ppz, 'String', valor_ppz);
       set(handles.txt_rmsz,'String',valor_rmsz);
       set(handles.txt_medz, 'String', acel_z);
       fzcap=fbz;
       if toc>=0 && toc<=1
           capz=[capz,fzcap];
           set(handles.txt_frecx, 'String', 'DATOS GUARDADOS');
       end
              else
               set(handles.txt_picoz, 'String', '');
               set(handles.txt_ppz, 'String', '');
               set(handles.txt_rmsz, 'String', '');
               set(handles.txt_medz, 'String', '');
       end
       axis([1 handles.xSamples -3 3]);
       drawnow;
       save('senalx',"capx");  
       save('senaly',"capy");
       save('senalz',"capz"); 
   t=8;
   Lx=numel(capx);
   Ly=numel(capy);
   Lz=numel(capz);
   Fsx=Lx/t;
   Fsy=Ly/t;
   Fsz=Lz/t;
   fx=Fsx/Lx*(0:Lx-1);
   fy=Fsy/Ly*(0:Ly-1);
   fz=Fsz/Lz*(0:Lz-1);
   ftx = fft(capx); 
   fty = fft(capy);
   ftz = fft(capz);

       axes(handles.axes2);
        
        if eje_x==1    
        plot(tiempo(1:contador-1),x(1:contador-1),Color='g'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
        else 
          
           
        end 
        if eje_y==1
        plot(tiempo(1:contador-1),y(1:contador-1), 'Color',[1, 0.5, 0]); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
       
        
      
        else
            
        end
        if eje_z==1
        plot(tiempo(1:contador-1),z(1:contador-1),'b'); hold on;
        xlabel(handles.axes2, 'TIEMPO');
        ylabel(handles.axes2, 'AMPLITUD');
      
        else
         
        end
       axis([1 handles.xSamples -3 3]);  
       pause(1);  
       end 
   if eje_x == 1 
figure;
plot(fx,abs(ftx),Color='g');
xlim([1 88]);
saveas(gcf, 'espectro_x.jpg');
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje X');  
end

if eje_y == 1
figure;
plot(fy, abs(fty), 'Color', 'r'); 
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Y');
saveas(gcf, 'espectro_y.jpg');
end

if eje_z == 1
figure;
plot(fz, abs(ftz), 'Color', 'b');  
xlim([1 88]);
xlabel('Frecuencia');
ylabel('FFT');
title('Espectro eje Z');
saveas(gcf, 'espectro_z.jpg');
end    
   end

end

function txt_fs_Callback(hObject, eventdata, handles)
handles.data1=get(hObject,'String');
handles.xSamples = str2double(handles.data1);

guidata(hObject,handles);


function txt_fs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt_medx_Callback(hObject, eventdata, handles)


function txt_medx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cmd_detener_Callback(hObject, eventdata, handles)
global stop;
stop = true;
cla(handles.axes1);
cla(handles.axes2);


function pm_sensor_Callback(hObject, eventdata, handles)
s1=get(hObject, 'String');

function pm_sensor_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_signal.
function pm_signal_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function pm_signal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
estado=get(hObject,'Value');
if estado==1
set(handles.pm_filtro, 'Enable', 'on');
set(handles.txt_fmax, 'Enable', 'on');
set(handles.txt_fmin, 'Enable', 'on');

else 
set(handles.pm_filtro, 'Enable', 'off');   
set(handles.txt_fmax, 'Enable', 'off');
set(handles.txt_fmin, 'Enable', 'off');
end
% --- Executes on selection change in pm_filtro.
function pm_filtro_Callback(hObject, eventdata, handles)
% hObject    handle to pm_filtro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_filtro contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_filtro


% --- Executes during object creation, after setting all properties.
function pm_filtro_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_fmax_Callback(hObject, eventdata, handles)

function txt_fmax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_fmin_Callback(hObject, eventdata, handles)

function txt_fmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_visual.
function pm_visual_Callback(hObject, eventdata, handles)
% hObject    handle to pm_visual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_visual contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_visual


% --- Executes during object creation, after setting all properties.
function pm_visual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_visual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_rmsx_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function txt_rmsx_CreateFcn(hObject, eventdata, handles)


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_picox_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function txt_picox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ppx_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ppx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ppx as text
%        str2double(get(hObject,'String')) returns contents of txt_ppx as a double


% --- Executes during object creation, after setting all properties.
function txt_ppx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ppx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_x.
function cb_x_Callback(hObject, eventdata, handles)
% ejes = get(handles.cb_x,'Value');
% if ejes == 1 
%     plot(y, 'LineWidth',0.5, Color='b');
% else
%     grid on;
% end


% --- Executes on button press in cb_y.
function cb_y_Callback(hObject, eventdata, handles)
% hObject    handle to cb_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_y


% --- Executes on button press in cb_z.
function cb_z_Callback(hObject, eventdata, handles)
% hObject    handle to cb_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_z



function txt_rmsz_Callback(hObject, eventdata, handles)
% hObject    handle to txt_rmsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_rmsz as text
%        str2double(get(hObject,'String')) returns contents of txt_rmsz as a double


% --- Executes during object creation, after setting all properties.
function txt_rmsz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_rmsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_picoz_Callback(hObject, eventdata, handles)
% hObject    handle to txt_picoz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_picoz as text
%        str2double(get(hObject,'String')) returns contents of txt_picoz as a double


% --- Executes during object creation, after setting all properties.
function txt_picoz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_picoz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ppz_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ppz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ppz as text
%        str2double(get(hObject,'String')) returns contents of txt_ppz as a double


% --- Executes during object creation, after setting all properties.
function txt_ppz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ppz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_rmsy_Callback(hObject, eventdata, handles)
% hObject    handle to txt_rmsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_rmsy as text
%        str2double(get(hObject,'String')) returns contents of txt_rmsy as a double



% --- Executes during object creation, after setting all properties.
function txt_rmsy_CreateFcn(hObject, ~, handles)
% hObject    handle to txt_rmsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function txt_picoy_Callback(hObject, ~, handles)
% hObject    handle to txt_picoy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_picoy as text
%        str2double(get(hObject,'String')) returns contents of txt_picoy as a double


% --- Executes during object creation, after setting all properties.
function txt_picoy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_picoy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ppy_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ppy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ppy as text
%        str2double(get(hObject,'String')) returns contents of txt_ppy as a double


% --- Executes during object creation, after setting all properties.
function txt_ppy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ppy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_magx_Callback(hObject, eventdata, handles)
% hObject    handle to txt_magx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_magx as text
%        str2double(get(hObject,'String')) returns contents of txt_magx as a double


% --- Executes during object creation, after setting all properties.
function txt_magx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_magx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_frecx_Callback(hObject, eventdata, handles)
% hObject    handle to txt_frecx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_frecx as text
%        str2double(get(hObject,'String')) returns contents of txt_frecx as a double


% --- Executes during object creation, after setting all properties.
function txt_frecx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_frecx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_medz_Callback(hObject, eventdata, handles)
% hObject    handle to txt_medz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_medz as text
%        str2double(get(hObject,'String')) returns contents of txt_medz as a double


% --- Executes during object creation, after setting all properties.
function txt_medz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_medz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_medy_Callback(hObject, eventdata, handles)
% hObject    handle to txt_medy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_medy as text
%        str2double(get(hObject,'String')) returns contents of txt_medy as a double


% --- Executes during object creation, after setting all properties.
function txt_medy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_medy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_magy_Callback(hObject, eventdata, handles)
% hObject    handle to txt_magy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_magy as text
%        str2double(get(hObject,'String')) returns contents of txt_magy as a double


% --- Executes during object creation, after setting all properties.
function txt_magy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_magy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_magz_Callback(hObject, eventdata, handles)
% hObject    handle to txt_magz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_magz as text
%        str2double(get(hObject,'String')) returns contents of txt_magz as a double



% --- Executes during object creation, after setting all properties.
function txt_magz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_magz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
