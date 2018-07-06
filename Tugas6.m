function varargout = Tugas6(varargin)
% TUGAS6 MATLAB code for Tugas6.fig
%      TUGAS6, by itself, creates a new TUGAS6 or raises the existing
%      singleton*.
%
%      H = TUGAS6 returns the handle to a new TUGAS6 or the handle to
%      the existing singleton*.
%
%      TUGAS6('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TUGAS6.M with the given input arguments.
%
%      TUGAS6('Property','Value',...) creates a new TUGAS6 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Tugas6_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Tugas6_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tugas6

% Last Modified by GUIDE v2.5 04-Apr-2018 11:38:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tugas6_OpeningFcn, ...
                   'gui_OutputFcn',  @Tugas6_OutputFcn, ...
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


% --- Executes just before Tugas6 is made visible.
function Tugas6_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Tugas6 (see VARARGIN)

% Choose default command line output for Tugas6
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tugas6 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tugas6_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnbrowse.
function btnbrowse_Callback(hObject, eventdata, handles)
% hObject    handle to btnbrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[namafile, formatfile] = uigetfile({'*.jpg','*.png'}, 'Pilih Gambar'); % memilih gambar
global image
image = imread([formatfile, namafile]); % membaca gambar
guidata(hObject, handles);
axes(handles.axes1); % memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); % memunculkan gambar


% --- Executes on button press in btnhistogram.
function btnhistogram_Callback(hObject, eventdata, handles)
% hObject    handle to btnhistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
R = image(:,:,1);
G = image(:,:,2);
B = image(:,:,3);

% red
Red = cat(3,R,G*0,B*0);
[rows,cols] = size(R); 

% histogram array
myhist = zeros(256,1);
for k = 0 : 255
    myhist(k+1) = numel(find(R == k)); % number of element dimana ada gray level = 'k'
end
figure, stem(myhist,'r');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on
 
% Green 
Green = cat(3,R*0,G,B*0);
[rows,cols] = size(G); 

% histogram array
myhist = zeros(256,1);
for k = 0 : 255
    myhist(k+1) = numel(find(G == k)); % number of element dimana ada gray level = 'k'
end
figure, stem(myhist,'g');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on
 
% blue 
Blue = cat(3,R*0,G*0,B);
[rows,cols] = size(B); 

% histogram array
myhist = zeros(256,1);
for k = 0 : 255
    myhist(k+1) = numel(find(B == k)); % number of element dimana ada gray level = 'k'
end
figure, stem(myhist,'b');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on

% --- Executes on button press in btnsembilan.
function btnsembilan_Callback(hObject, eventdata, handles)
% hObject    handle to btnsembilan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[rowsi,colsi,z]= size(image); 

sudut = 90;
rads = 2 * pi * sudut / 360;  

% menghitung dimesi array agar gambar yang diputar akan sesuai dengan yg didalamnya
% menggunakan absolute agar bisa mendapatkan nilai positif 
rowsf = ceil(rowsi * abs(cos(rads)) + colsi * abs(sin(rads)));                      
colsf = ceil(rowsi * abs(sin(rads)) + colsi * abs(cos(rads)));                     

% menentukan array dengan parameter yang dihitung dan isi array dengan angka nol
C = uint8(zeros([rowsf colsf 3 ]));

% menghitung center of original and final image
xo = ceil(rowsi / 2);                                                            
yo = ceil(colsi / 2);

midx = ceil((size(C,1)) / 2);
midy = ceil((size(C,2)) / 2);

% menghitung koordinat yang sesuai dari pixel A
% untuk setiap pixel C dan intensitasnya akan diberikan setelah diperiksa
for i = 1 : size(C,1)
    for j = 1 : size(C,2)                                                       
         x = (i - midx) * cos(rads) + (j - midy) * sin(rads);                                       
         y = -(i - midx) * sin(rads) +(j - midy) * cos(rads);                             
         x = round(x) + xo;
         y = round(y) + yo;
         if (x >= 1 && y >= 1 && x <= size(image,1) &&  y <= size(image,2)) 
              C(i,j,:) = image(x,y,:);  
         end
    end
end
image = C;
axes(handles.axes1); % memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); % memunculkan gambar


% --- Executes on button press in btndelapan.
function btndelapan_Callback(hObject, eventdata, handles)
% hObject    handle to btndelapan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[rowsi,colsi,z]= size(image); 

sudut = 180;
rads = 2 * pi * sudut / 360;  

% menghitung dimesi array agar gambar yang diputar akan sesuai dengan yg didalamnya
% menggunakan absolute agar bisa mendapatkan nilai positif 
rowsf = ceil(rowsi * abs(cos(rads)) + colsi * abs(sin(rads)));                      
colsf = ceil(rowsi * abs(sin(rads)) + colsi * abs(cos(rads)));                     

% menentukan array dengan parameter yang dihitung dan isi array dengan angka nol
C = uint8(zeros([rowsf colsf 3 ]));

% menghitung center of original and final image
xo = ceil(rowsi / 2);                                                            
yo = ceil(colsi / 2);

midx = ceil((size(C,1)) / 2);
midy = ceil((size(C,2)) / 2);

% menghitung koordinat yang sesuai dari pixel A
% untuk setiap pixel C dan intensitasnya akan diberikan setelah diperiksa
for i = 1 : size(C,1)
    for j = 1 : size(C,2)                                                       
         x = (i - midx) * cos(rads) + (j - midy) * sin(rads);                                       
         y = -(i - midx) * sin(rads) +(j - midy) * cos(rads);                             
         x = round(x) + xo;
         y = round(y) + yo;
         if (x >= 1 && y >= 1 && x <= size(image,1) &&  y <= size(image,2)) 
              C(i,j,:) = image(x,y,:);  
         end
    end
end
image = C;
axes(handles.axes1); % memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); % memunculkan gambar


% --- Executes on button press in btntujuh.
function btntujuh_Callback(hObject, eventdata, handles)
% hObject    handle to btntujuh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[rowsi,colsi,z]= size(image); 

sudut = 270;
rads = 2 * pi * sudut / 360;  

% menghitung dimesi array agar gambar yang diputar akan sesuai dengan yg didalamnya
% menggunakan absolute agar bisa mendapatkan nilai positif 
rowsf = ceil(rowsi * abs(cos(rads)) + colsi * abs(sin(rads)));                      
colsf = ceil(rowsi * abs(sin(rads)) + colsi * abs(cos(rads)));                     

% menentukan array dengan parameter yang dihitung dan isi array dengan angka nol
C = uint8(zeros([rowsf colsf 3 ]));

% menghitung center of original and final image
xo = ceil(rowsi / 2);                                                            
yo = ceil(colsi / 2);

midx = ceil((size(C,1)) / 2);
midy = ceil((size(C,2)) / 2);

% menghitung koordinat yang sesuai dari pixel A
% untuk setiap pixel C dan intensitasnya akan diberikan setelah diperiksa
for i = 1 : size(C,1)
    for j = 1 : size(C,2)                                                       
         x = (i - midx) * cos(rads) + (j - midy) * sin(rads);                                       
         y = -(i - midx) * sin(rads) +(j - midy) * cos(rads);                             
         x = round(x) + xo;
         y = round(y) + yo;
         if (x >= 1 && y >= 1 && x <= size(image,1) &&  y <= size(image,2)) 
              C(i,j,:) = image(x,y,:);  
         end
    end
end
image = C;
axes(handles.axes1); % memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); % memunculkan gambar


% --- Executes on button press in btngray.
function btngray_Callback(hObject, eventdata, handles)
% hObject    handle to btngray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image
merah =image(:,:,1); % hanya berisi piksel warna merah
ijo = image(:,:,2);% hanya berisi piksel warna merah
biru = image(:,:,3); % hanya berisi piksel warna merah

%ubah warna ke abu2
abu = (0.3 * merah) + (0.5 * ijo) + (0.2 * biru);
image = abu;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnzoomin.
function btnzoomin_Callback(hObject, eventdata, handles)
% hObject    handle to btnzoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s = size(image);
f1 = 200;
s2=s*f1;
k=1;
l=1;
for (i=1:f1:s2)
    for( j=1:f1:s2)
        C(i,j)= image(k,l);
        l=l+1;
    end
    l=1;
    k=k+1;
end

for (i=1:f1:s2)
    for (j=2:f1:s2-1)
        C(i,j)= [C(i,j-1)+ C(i, j+1)]*0.5;
    end
end

for(j=1:f1:s2)
    for(i=2:f1:s2-1)
        C(i,j)=[C(i-1,j)+C(i+1,j)]*0.5;
    end
end

for (i=2:f1:s2-1)
    for (j=2:f1:s2-1)
        C(i,j)= [C(i,j-1)+ C(i, j+1)]*0.5;
    end
end
image = C;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnzoomout.
function btnzoomout_Callback(hObject, eventdata, handles)
% hObject    handle to btnzoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
% ukuran gambar 512X512 --> 256x256 piksel
area=uint8(zeros(256,256));
% menduplikasikan data citra asli
d=zeros(512,512);
for i=1:512;
    for j=1:512;
        d(i,j)=image(i,j);
    end
end
%kompresi dari 512x512 ke 256x256 piksel
for b_asli=1:256
    for k_asli=1:256
        temp=0; dummy=0;
        for b_baru=1:2
            for k_baru=1:2
                dummyb=((b_asli-1)*2 + b_baru);
                dummyk= ((k_asli-1)*2 + k_baru);
                dummy=d(dummyb,dummyk);
                temp=temp+dummy;
            end
        end
        temp=round(temp/4);
        area(b_asli,k_asli)=temp;
    end
end
image = area;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btninverse.
function btninverse_Callback(hObject, eventdata, handles)
% hObject    handle to btninverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
nilai = uint8(image);
baris = size(nilai,1);
kolom = size(nilai,2);
for i = 1 : baris
    for j = 1 : kolom
        temp = 255 - nilai(i,j);
        nilai(i,j) = temp;
    end
end
image = nilai;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btntambah.
function btntambah_Callback(hObject, eventdata, handles)
% hObject    handle to btntambah (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
nilai = uint8(image);
baris = size(nilai,1);
kolom = size(nilai,2);
for i = 1 : baris
    for j = 1 : kolom
        temp = nilai(i,j) + 20;
        nilai(i,j) = temp;
    end
end
image = nilai;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnkali.
function btnkali_Callback(hObject, eventdata, handles)
% hObject    handle to btnkali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
nilai = uint8(image);
baris = size(nilai,1);
kolom = size(nilai,2);
for i = 1 : baris
    for j = 1 : kolom
        temp = nilai(i,j) * 2;
        nilai(i,j) = temp;
    end
end
image = nilai;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnkurang.
function btnkurang_Callback(hObject, eventdata, handles)
% hObject    handle to btnkurang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
nilai = uint8(image);
baris = size(nilai,1);
kolom = size(nilai,2);
for i = 1 : baris
    for j = 1 : kolom
        temp = nilai(i,j) - 10;
        nilai(i,j) = temp;
    end
end
image = nilai;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnbagi.
function btnbagi_Callback(hObject, eventdata, handles)
% hObject    handle to btnbagi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
nilai = uint8(image);
baris = size(nilai,1);
kolom = size(nilai,2);
for i = 1 : baris
    for j = 1 : kolom
        temp = nilai(i,j) / 2;
        nilai(i,j) = temp;
    end
end
image = nilai;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnhorizontal.
function btnhorizontal_Callback(hObject, eventdata, handles)
% hObject    handle to btnhorizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
ind = 1; %inisialisasi index
for i = size(image,2) : -1 : 1
    new(:, ind, :) = image(:, i, :);
    ind = ind + 1;
end
image = new; %flipdim(image, 2);
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnvertikal.
function btnvertikal_Callback(hObject, eventdata, handles)
% hObject    handle to btnvertikal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
ind = 1; %inisialisasi index
for i = size(image,1) : -1 : 1
    new(ind, :, :) = image(i, :, :);
    ind = ind + 1;
end
image = new; %flipdim(image, 1);
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btncrop.
function btncrop_Callback(hObject, eventdata, handles)
% hObject    handle to btncrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
x_TopLeft = 64;
y_TopLeft = 64;
x_BottomRight = 192; 
y_BottomRight = 192;

M = x_BottomRight-x_TopLeft+1;
N = y_BottomRight-y_TopLeft+1; 

for i = 1 : M
    for j = 1 : N
        B(i, j, :) = image(x_TopLeft+i, y_TopLeft+j, :);
    end
end
image = B;
axes(handles.axes1); %memilih plotori sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --- Executes on button press in btnhist.
function btnhist_Callback(hObject, eventdata, handles)
% hObject    handle to btnhist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
R = image(:,:,1);
G = image(:,:,2);
B = image(:,:,3);

% red
Red = cat(3,R,G*0,B*0);
[rows,cols] = size(R); 

% histogram array
myhist = zeros(256,1);
for k = 0 : 255
    myhist(k+1) = numel(find(R == k)); % number of element dimana ada gray level = 'k'
end
figure, stem(myhist,'r');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on
 
% Green 
Green = cat(3,R*0,G,B*0);
[rows,cols] = size(G); 

% histogram array
myhist = zeros(256,1);
for k = 0 : 255
    myhist(k+1) = numel(find(G == k)); % number of element dimana ada gray level = 'k'
end
figure, stem(myhist,'g');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on
 
% blue 
Blue = cat(3,R*0,G*0,B);
[rows,cols] = size(B); 

% histogram array
myhist = zeros(256,1);
for k = 0 : 255
    myhist(k+1) = numel(find(B == k)); % number of element dimana ada gray level = 'k'
end
figure, stem(myhist,'b');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on


% --- Executes on button press in btnsembilan.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to btnsembilan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btndelapan.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to btndelapan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btntujuh.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to btntujuh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnedge.
function btnedge_Callback(hObject, eventdata, handles)
% hObject    handle to btnedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[row, coll, z] = size(image);
gambarbaru = double(image);
y = [1 1 1;
    1 -8 1;
    1 1 1];
for i = 2 : row - 2
    for j = 2 : coll - 2
        hasil = gambarbaru(i-1,j-1)*y(1,1)+gambarbaru(i,j-1)*y(2,1)...
            + gambarbaru(i+1,j-1)*y(3,1)+gambarbaru(i-1,j)*y(1,2)...
            + gambarbaru(i,j)*y(2,2)+gambarbaru(i+1,j)*y(3,2)...
            + gambarbaru(i-1,j+1)*y(1,3)+gambarbaru(i,j+1)*y(2,3)...
            + gambarbaru(i+1,j+1)*y(3,3);
        gmb(i-1,j-1) = hasil;
    end
end
figure
imshow(uint8(gmb));


% --- Executes on button press in btnblur.
function btnblur_Callback(hObject, eventdata, handles)
% hObject    handle to btnblur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[row, coll, z] = size(image);
gambarbaru = double(image);
y = [1/9 1/9 1/9;
    1/9 1/9 1/9;
    1/9 1/9 1/9];
for i = 2 : row - 2
    for j = 2 : coll - 2
        hasil = gambarbaru(i-1,j-1,:)*y(1,1,:)+gambarbaru(i,j-1,:)*y(2,1,:)...
            + gambarbaru(i+1,j-1,:)*y(3,1,:)+gambarbaru(i-1,j,:)*y(1,2,:)...
            + gambarbaru(i,j,:)*y(2,2,:)+gambarbaru(i+1,j,:)*y(3,2,:)...
            + gambarbaru(i-1,j+1,:)*y(1,3,:)+gambarbaru(i,j+1,:)*y(2,3,:)...
            + gambarbaru(i+1,j+1,:)*y(3,3,:);
        gmb(i-1,j-1,:) = hasil;
    end
end
figure
imshow(uint8(gmb));


% --- Executes on button press in btnsharp.
function btnsharp_Callback(hObject, eventdata, handles)
% hObject    handle to btnsharp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[row, coll, z] = size(image);
gambarbaru = double(image);
y = [0 -1 0;
    -1 5 -1;
    0 -1 0];
for i = 2 : row -2
    for j = 2 : coll -2
        hasil = gambarbaru(i-1,j-1,:)*y(1,1,:)+gambarbaru(i,j-1,:)*y(2,1,:)...
            + gambarbaru(i+1,j-1,:)*y(3,1,:)+gambarbaru(i-1,j,:)*y(1,2,:)...
            + gambarbaru(i,j,:)*y(2,2,:)+gambarbaru(i+1,j,:)*y(3,2,:)...
            + gambarbaru(i-1,j+1,:)*y(1,3,:)+gambarbaru(i,j+1,:)*y(2,3,:)...
            + gambarbaru(i+1,j+1,:)*y(3,3,:);
        gmb(i-1,j-1,:) = hasil;
    end
end
figure
imshow(uint8(gmb));



% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnnmean.
function btnnmean_Callback(hObject, eventdata, handles)
% hObject    handle to btnnmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[row, col, z] = size (image);
gambarbaru = double(image);
for i = 1 : row %h
    for j = 1 : col %w
        if (i == row) 
            if (i == row && j == col)
                hasil = (gambarbaru(i,j,:) + 0 + 0 + 0)/4;
            else
                hasil = (gambarbaru(i,j,:) + gambarbaru(i,j+1,:) + 0 + 0)/4;
            end
        else if (j == col)
            if (j == col && i == row)
                hasil = (gambarbaru(i,j,:) + 0 + 0 + 0)/4;
            else
                hasil = (gambarbaru(i,j,:) + 0 + gambarbaru(i+1,j,:) + 0)/4; 
            end
        else
            hasil = (gambarbaru(i,j,:) + gambarbaru(i,j+1,:) + gambarbaru(i+1,j,:) + gambarbaru(i+1,j+1,:))/4;
            end 
        end
        gmb(i,j,:) = hasil;
    end
end
figure
imshow(uint8(gmb));


% --- Executes on button press in btnnmedian.
function btnnmedian_Callback(hObject, eventdata, handles)
% hObject    handle to btnnmedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[row, col, z] = size (image);
gambarbaru = double(image);
for i = 1 : row %h
    for j = 1 : col %w
        if (i == row) 
            if (i == row && j == col)
                hasil = 0;
            else
                a = gambarbaru(i,j,:);b = gambarbaru(i,j+1,:);
                g = [a b]; s = sort(g); n = s(1,1);
                hasil = median([0 n]);
            end
        else if (j == col)
            if (j == col && i == row)
                hasil = median([gambarbaru(i,j,:), 0, 0, 0]);
            else
                a = gambarbaru(i,j,:);b = gambarbaru(i+1,j,:);
                g = [a b]; s = sort(g); n = s(1,1);
                hasil = median([0 n]);
            end
        else
            hasil = median([gambarbaru(i,j,:), gambarbaru(i,j+1,:), gambarbaru(i+1,j,:), gambarbaru(i+1,j+1,:)]);
            end 
        end
        gmb(i,j,:) = hasil;
    end
end
figure
imshow(uint8(gmb));



% --- Executes on button press in btnnmodus.
function btnnmodus_Callback(hObject, eventdata, handles)
% hObject    handle to btnnmodus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[row, col, z] = size (image);
gambarbaru = double(image);
for i = 1 : row %h
    for j = 1 : col %w
        if (i == row) 
            if (i == row && j == col)
                hasil = 0;
            else
                hasil = 0;
            end
        else if (j == col)
            if (j == col && i == row)
                hasil = mode([gambarbaru(i,j,:), 0, 0, 0]);
            else
                hasil = 0;
            end
        else
            hasil = mode([gambarbaru(i,j,:), gambarbaru(i,j+1,:), gambarbaru(i+1,j,:), gambarbaru(i+1,j+1,:)]);
            end 
        end
        gmb(i,j,:) = hasil;
    end
end
figure
imshow(uint8(gmb));
