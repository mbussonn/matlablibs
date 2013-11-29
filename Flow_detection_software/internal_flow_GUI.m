function varargout = internal_flow_GUI(varargin)
% INTERNAL_FLOW_GUI M-file for internal_flow_GUI.fig
%      INTERNAL_FLOW_GUI, by itself, creates a new INTERNAL_FLOW_GUI or raises the existing
%      singleton*.
%
%      H = INTERNAL_FLOW_GUI returns the handle to a new INTERNAL_FLOW_GUI or the handle to
%      the existing singleton*.
%
%      INTERNAL_FLOW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERNAL_FLOW_GUI.M with the given input arguments.
%
%      INTERNAL_FLOW_GUI('Property','Value',...) creates a new INTERNAL_FLOW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before internal_flow_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to internal_flow_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help internal_flow_GUI

% Last Modified by GUIDE v2.5 23-Dec-2011 17:38:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @internal_flow_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @internal_flow_GUI_OutputFcn, ...
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


% --- Executes just before internal_flow_GUI is made visible.
function internal_flow_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to internal_flow_GUI (see VARARGIN)

% Choose default command line output for internal_flow_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes internal_flow_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = internal_flow_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
files=handles.files;
d_path=handles.d_path;
    i=round(get(handles.slider1,'Value'));
    im=imread([d_path,'\',files(i).name]);
     im=im(:,:,get(handles.listbox_color,'Value'));
    handles.im=im;      
    im2=imread([d_path,'\',files(i+1).name]);
     im2=im2(:,:,get(handles.listbox_color,'Value'));
    handles.im2=im2;
    
handles=detect_edge(hObject, eventdata, handles)
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function load_im_Callback(hObject, eventdata, handles)
% hObject    handle to load_im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider1,'Value',1)
d_path='E:\Science\data\benoit\7511_up_0920-1005'
d_path=uigetdir(d_path)
    files=dir([d_path,'\*.png']);
    curr=round(get(handles.slider1,'Value'));
    i=curr;
    im=imread([d_path,'\',files(i).name]);
    handles.im=im(:,:,get(handles.listbox_color,'Value'));      
    im2=imread([d_path,'\',files(i+1).name]);
    handles.im2=im2(:,:,get(handles.listbox_color,'Value'));

    set(handles.slider1,'Max',length(files)-1);
    slider_step(1) = 1/(length(files)-1);
    slider_step(2) = 1/(length(files)-1);
    set(handles.slider1,'sliderstep',slider_step);
    


axes(handles.image);
imshow(handles.im);
handles.curr=curr;
handles.d_path=d_path;
handles.files=files;
handles=detect_edge(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=correlate_im(hObject, eventdata, handles);
guidata(hObject, handles);


function s_source_Callback(hObject, eventdata, handles)
% hObject    handle to s_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_source as text
%        str2double(get(hObject,'String')) returns contents of s_source as a double


% --- Executes during object creation, after setting all properties.
function s_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function s_search_Callback(hObject, eventdata, handles)
% hObject    handle to s_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_search as text
%        str2double(get(hObject,'String')) returns contents of s_search as a double


% --- Executes during object creation, after setting all properties.
function s_search_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g_dist_Callback(hObject, eventdata, handles)
% hObject    handle to g_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g_dist as text
%        str2double(get(hObject,'String')) returns contents of g_dist as a double


% --- Executes during object creation, after setting all properties.
function g_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to c_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_thresh as text
%        str2double(get(hObject,'String')) returns contents of c_thresh as a double


% --- Executes during object creation, after setting all properties.
function c_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flow_struct=save_actual_im(hObject, eventdata, handles)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    files=handles.files;
    d_path=handles.d_path;

for i=1:length(files)-1
    set(handles.slider1,'value',i)
    im=imread([d_path,'\',files(i).name]);
    im=im(:,:,get(handles.listbox_color,'Value'));
    handles.im=im;      
    im2=imread([d_path,'\',files(i+1).name]);
    im2=im2(:,:,get(handles.listbox_color,'Value'));
    handles.im2=im2;

    handles.curr=i;
    guidata(hObject, handles);
    handles=detect_edge(hObject, eventdata, handles)
    handles=correlate_im(hObject, eventdata, handles)
    guidata(hObject, handles);
    flow_struct=save_actual_im(hObject, eventdata, handles)
end


function size_er_Callback(hObject, eventdata, handles)
% hObject    handle to size_er (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_er as text
%        str2double(get(hObject,'String')) returns contents of size_er as a double
handles=detect_edge(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function size_er_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_er (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function size_dil_Callback(hObject, eventdata, handles)
% hObject    handle to size_dil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_dil as text
%        str2double(get(hObject,'String')) returns contents of size_dil as a double

handles=detect_edge(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function size_dil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_dil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles=detect_edge(hObject, eventdata, handles)
guidata(hObject, handles);
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



function remove_edge_p_Callback(hObject, eventdata, handles)
% hObject    handle to remove_edge_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of remove_edge_p as text
%        str2double(get(hObject,'String')) returns contents of remove_edge_p as a double

handles=detect_edge(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function remove_edge_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to remove_edge_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles=correlate_im(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gd=str2num(get(handles.g_dist,'String'));
gs=str2num(get(handles.s_source,'String'));
ss=str2num(get(handles.s_search,'String'));
c_th=str2num(get(handles.c_thresh,'String'));


curr=round(get(handles.slider1,'Value'));
axes(handles.image);
imshow(uint8(handles.im_e))
handles.curr;

im=handles.im;

im2=handles.im2;
BW2=handles.BW2;
se = strel('disk',gs);
BW2=imerode(BW2,se);
n=1;
for i=1:gd:size(im,1)
    for j=1:gd:size(im,2)
        if BW2(i,j)==1
            try 
                sub_roi=im(i-gs:i+gs,j-gs:j+gs);
            catch
                continue
            end
            try
                sub_area=im2(i-ss:i+ss,j-ss:j+ss);
            catch
                continue
            end
            try
                c = normxcorr2(sub_roi,sub_area);
            catch
                continue
            end
                [xroi,yroi]=size(sub_roi);
            [xarea,yarea]=size(sub_area);
            c=c(xroi:xarea,yroi:yarea);
            % figure, surf(c), shading flat
            % offset found by correlation
            [max_c, imax] = max(c(:));
            [ypeak, xpeak] = ind2sub(size(c),imax(1));
            [xc,yc]=size(c);
%             x_p(n)=(j);
%             y_p(n)=(i);
%             x_v(n)=xpeak-(xc+1)/2-1;
%             y_v(n)=ypeak-(yc+1)/2-1;
            x_vector(i,j)=xpeak-(xc+1)/2;
            y_vector(i,j)=ypeak-(yc+1)/2;
            corr_val(i,j)=max_c;
            n=n+1;
            
            if max_c>=c_th
                x_p(n)=(j);
                y_p(n)=(i);
                x_v(n)=xpeak-(xc+1)/2;
                y_v(n)=ypeak-(yc+1)/2;
                c_val(n)=max_c;
%                 axes(handles.gran);
%                 surf(c)
%                 pause(.1)
                %max_c
            end
            
        end
    end
end
hold on
quiver(x_p,y_p,x_v,y_v,'r');

hold off

handles.x_pos=x_p;
handles.y_pos=y_p;
handles.x_vec=x_v;
handles.y_vec=y_v;
handles.c_val=c_val;

guidata(hObject, handles);
%now I create the grid


% --- Executes on slider movement.
function handles=detect_edge(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


curr=round(get(handles.slider1,'Value'));
set(handles.text1,'String',['Image: ',num2str(curr)]);
axes(handles.image);
handles.curr;
im=handles.im;
im=im;
imshow(im);
guidata(hObject, handles);

axes(handles.edge);
level = graythresh(im);

e_meth_l=get(handles.popupmenu1,'String');
e_meth_v=get(handles.popupmenu1,'Value');
e_meth=e_meth_l{e_meth_v};


BW1 = edge(im,e_meth);

imshow(BW1);

axes(handles.gran);

[L,num] = bwlabel(BW1);
stats = regionprops(L,'Area');
allArea = [stats.Area];
%now_we remove small particles
BW2=BW1;
size_cut=str2num(get(handles.remove_edge_p,'String'));
li=find(allArea<size_cut);
for i=1:length(li)
    BW2(find(L==li(i)))=0;
end

se_size=str2num(get(handles.size_dil,'String'));
se2_size=str2num(get(handles.size_er,'String'));
se = strel('disk',se_size);
se2 = strel('disk',se2_size);
BW2=imdilate(BW2,se);
BW2f=imfill(BW2,4,'holes');
BW2=imerode(BW2f,se2);
%imshow(BW2);
imshow(BW2);

axes(handles.image);
BW2_i=edge(double(BW2));
BW2_i=-1*BW2_i+1;
im_e=double(im).*double(BW2_i);
imshow(uint8(im_e));

handles.BW2=BW2;
handles.im_e=im_e;
guidata(hObject, handles);


function flow_struct=save_actual_im(hObject, eventdata, handles)

im=handles.im;
files=handles.files;
curr=round(get(handles.slider1,'value'))
%First I arrange the structure
flow_struct.d_path=handles.d_path;
flow_struct.curr=curr;
flow_struct.im=im;
flow_struct.x_pos=handles.x_pos;
flow_struct.y_pos=handles.y_pos;
flow_struct.x_vec=handles.x_vec;
flow_struct.y_vec=handles.y_vec;
flow_struct.c_val=handles.c_val;
flow_struct.file=files(curr).name;
flow_struct.edge_im=handles.BW2;

e_meth_l=get(handles.popupmenu1,'String')
e_meth_v=get(handles.popupmenu1,'Value')
e_meth=e_meth_l{e_meth_v}
flow_struct.edge_meth=e_meth;

flow_struct.size_cut=str2num(get(handles.remove_edge_p,'String'));
flow_struct.se_size=str2num(get(handles.size_dil,'String'));
flow_struct.se2_size=str2num(get(handles.size_er,'String'));


flow_struct.gd=str2num(get(handles.g_dist,'String'));
flow_struct.gs=str2num(get(handles.s_source,'String'));
flow_struct.ss=str2num(get(handles.s_search,'String'));
flow_struct.c_th=str2num(get(handles.c_thresh,'String'));

[status,message,messageid] = mkdir(handles.d_path,'data');
[status,message,messageid] = mkdir(handles.d_path,'image');
%this is not done
fil=files(curr).name;
axes(handles.image);
save([handles.d_path,'\data\',fil(1:end-3),'mat'],'flow_struct')
%print([handles.d_path,'\image\',fil(1:end-4),'_results.png'],'-fhandles.image','-dpng','-r300')
%saveas(handles.image,[handles.d_path,'\image\',fil(1:end-4),'_results.png'])
newfigure = figure;
imshow(uint8(handles.im_e))
hold on
quiver(handles.x_pos,handles.y_pos,handles.x_vec,handles.y_vec,'r');
% print(gcf,[handles.d_path,'\image\',fil(1:end-4),'_results.png'],'-dpng','-r150')
saveas(gcf,[handles.d_path,'\image\',fil(1:end-4),'_results.png'])
% close the new figure
close(newfigure) 


% --- Executes on selection change in listbox_color.
function listbox_color_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_color


% --- Executes during object creation, after setting all properties.
function listbox_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kernel_ee_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_ee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_ee as text
%        str2double(get(hObject,'String')) returns contents of kernel_ee as a double


% --- Executes during object creation, after setting all properties.
function kernel_ee_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_ee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kernel_size_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_size as text
%        str2double(get(hObject,'String')) returns contents of kernel_size as a double


% --- Executes during object creation, after setting all properties.
function kernel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=interpolate_flow(hObject,handles)
%Here we will apply an gaussian interpolation. Each measured point will
%have a gaussian field associated with it. The field will have the decay
%and the size as defined in the GUI. In the end, the average of each values
%is calculated and the image field is displayed and stored in the handle.
%Later it will be also put into the save structure that is stored if
%desired.
