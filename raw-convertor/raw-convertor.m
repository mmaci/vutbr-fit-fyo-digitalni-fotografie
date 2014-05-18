function varargout = raw(varargin)
% RAW MATLAB code for raw.fig
%      RAW, by itself, creates a new RAW or raises the existing
%      singleton*.
%
%      H = RAW returns the handle to a new RAW or the handle to
%      the existing singleton*.
%
%      RAW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAW.M with the given input arguments.
%
%      RAW('Property','Value',...) creates a new RAW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before raw_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to raw_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help raw

% Last Modified by GUIDE v2.5 25-Mar-2014 01:42:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @raw_OpeningFcn, ...
                   'gui_OutputFcn',  @raw_OutputFcn, ...
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


% Applies (inverse) sRGB gamma correction to an RGB image.
function out = srgb_gamma(in)

out = zeros(size(in));
nl=in>0.0031308;
out(nl)=1.055*(in(nl).^(1/2.4)) - 0.055;
out(~nl)=12.92*in(~nl);


% Applies a color space conversion matrix
function out = apply_cspace_matrix(im, cmatrix)

r = cmatrix(1,1)*im(:,:,1)+cmatrix(1,2)*im(:,:,2)+cmatrix(1,3)*im(:,:,3);
g = cmatrix(2,1)*im(:,:,1)+cmatrix(2,2)*im(:,:,2)+cmatrix(2,3)*im(:,:,3);
b = cmatrix(3,1)*im(:,:,1)+cmatrix(3,2)*im(:,:,2)+cmatrix(3,3)*im(:,:,3);

out = cat(3,r,g,b);


% --- Executes just before raw is made visible.
function raw_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to raw (see VARARGIN)




% Choose default command line output for raw
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes raw wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = raw_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes when selected object is changed in uipanel2.
function bgroup_SelectionChangeFcn(hObject, eventdata, handles)
global meta cfa cropped balanced wb demosaiced srgb postprocess stage;

set(handles.slider_saturation, 'Enable' ,'off');
set(handles.slider_contrast, 'Enable' ,'off');
set(handles.slider_lightness, 'Enable' ,'off');

set(handles.b_he, 'Enable' ,'off');
set(handles.b_restore, 'Enable' ,'off');
set(handles.b_gamma, 'Enable' ,'off');
set(handles.b_lrotate, 'Enable' ,'off');
set(handles.b_rrotate, 'Enable' ,'off');  

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    % initial stage
    case 'rb_cfa'
        postprocess = cfa;
        axes(handles.display)
        imshow(cfa);
        
    % crops image to active areas
    case 'rb_cropped'
        if stage < 1
            x_0 = meta.SubIFDs{1}.ActiveArea(2)+1;
            width = meta.SubIFDs{1}.DefaultCropSize(1);
            y_0 = meta.SubIFDs{1}.ActiveArea(1)+1;
            height = meta.SubIFDs{1}.DefaultCropSize(2);

            cropped = cfa(y_0:y_0+height,x_0:x_0+width);                   
            
            stage = 1;            
        end
        set(handles.rb_balanced, 'Enable' ,'on');
        postprocess = cropped;
        axes(handles.display);
        imshow(cropped);        
        
    % balances black and white levels    
    case 'rb_balanced'
        if stage < 2
            temp = double(cropped);
            black = meta.SubIFDs{1}.BlackLevel;
            white = meta.SubIFDs{1}.WhiteLevel;

            temp = temp - black;
            temp = temp ./ (white-black);

            balanced = max(0, min(temp, 1));
            stage = 2;
            
        end
        set(handles.rb_wb, 'Enable' ,'on');
        postprocess = balanced;
        axes(handles.display);
        imshow(balanced);        
        
    % white balance
    case 'rb_wb'  
        if stage < 3
            wb_coeffs = meta.AsShotNeutral .^ (-1);
            wb_coeffs = wb_coeffs ./ wb_coeffs(2);

            % note that CFA must be:
            % R G
            % G B
            wb_mask = ones(size(balanced));
            % every odd element both x-wise, y-wise
            wb_mask(1:2:end, 1:2:end) = wb_coeffs(1);
            % every even element both x-wise, y-wise
            wb_mask(2:2:end, 2:2:end) = wb_coeffs(3); 

            wb = balanced .* wb_mask;
            
            stage = 3;
            
        end
        set(handles.rb_demosaiced, 'Enable' ,'on');
        postprocess = wb;
        axes(handles.display);
        imshow(wb);                
        
    % demosaicing
    % must be applied to rggb sensors
    case 'rb_demosaiced'
        if stage < 4
            temp = uint16(wb .* (2^16));
            temp = demosaic(temp,'rggb');
            demosaiced = double(temp) ./ (2^16);
            
            stage = 4;
        end
        set(handles.rb_srgb, 'Enable' ,'on');
        postprocess = demosaiced;
        axes(handles.display);
        imshow(demosaiced);

    % conversion to srgb color space
    case 'rb_srgb'
        if stage < 5
            xyz_to_camera = meta.ColorMatrix1;
            xyz_to_camera = reshape(xyz_to_camera, 3, 3);

            % http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
            srgb_to_xyz = [0.4124564 0.3575761 0.1804375;
                           0.2126729 0.7151522 0.0721750;
                           0.0193339 0.1191920 0.9503041];

            srgb_to_camera = srgb_to_xyz * xyz_to_camera;
            srgb_to_camera = srgb_to_camera ./ repmat(sum(srgb_to_camera, 2), 1, 3);
            camera_to_srgb = srgb_to_camera ^(-1);

            temp = apply_cspace_matrix(demosaiced, camera_to_srgb);
            srgb = max(0,min(temp,1));
            stage = 5;                                    
        end
        
        % and now we are ready to process the image
        % thus we enable all the buttons
        set(handles.slider_saturation, 'value', 0.5);
        set(handles.slider_contrast, 'value', 0.5);
        set(handles.slider_lightness, 'value', 0.5);

        set(handles.slider_saturation, 'Enable' ,'on');
        set(handles.slider_contrast, 'Enable' ,'on');
        set(handles.slider_lightness, 'Enable' ,'on');                

        set(handles.b_he, 'Enable' ,'on');
        set(handles.b_restore, 'Enable' ,'on');
        set(handles.b_gamma, 'Enable' ,'on');
        set(handles.b_lrotate, 'Enable' ,'on');
        set(handles.b_rrotate, 'Enable' ,'on');
        
        postprocess = srgb;
        axes(handles.display);
        imshow(postprocess);        
 
    otherwise
        % Code for when there is no match.
end


% --- Executes on button press in loadFileButton.
% --- Loads a dng file
function loadFileButton_Callback(hObject, eventdata, handles)
    global meta cfa stage;

    set(handles.b_save, 'Enable' ,'on');
    
    set(handles.rb_balanced, 'Enable' ,'off');
    set(handles.rb_demosaiced, 'Enable' ,'off');
    set(handles.rb_wb, 'Enable' ,'off');
    set(handles.rb_srgb, 'Enable' ,'off');    
    set(handles.b_he, 'Enable' ,'off');
    set(handles.b_restore, 'Enable' ,'off');
    set(handles.b_gamma, 'Enable' ,'off');
    set(handles.b_lrotate, 'Enable' ,'off');
    set(handles.b_rrotate, 'Enable' ,'off');    


    [filename,filepath]=uigetfile({'*.*','All Files'},...
      'Load DNG File');
      cd(filepath);

    tiff = Tiff(filename, 'r');
    offsets = getTag(tiff, 'SubIFD');
    setSubDirectory(tiff, offsets(1));

    cfa = read(tiff);    
    close(tiff);
    
    meta = imfinfo(filename);
    set(handles.display, 'Visible' ,'on');
    axes(handles.display)
    imshow(cfa)

    set(handles.rb_cfa, 'Enable' ,'on');
    set(handles.rb_cropped, 'Enable' ,'on');

    stage = 0;


% --- Executes on slider movement.
% --- Saturates colors
function slider_saturation_Callback(hObject, eventdata, handles)
    global postprocess;

    value = get(hObject, 'value');

    postprocess = rgb2hsv(postprocess);
    multiplier = 2 * value;
    postprocess(:, :, 2) = min(1, max(0, postprocess(:, :, 2) * multiplier));
    postprocess = hsv2rgb(postprocess);
    axes(handles.display);    
    imshow(postprocess);


% --- Executes during object creation, after setting all properties.
function slider_saturation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_saturation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
% --- Enhances or decreases contrast
function slider_contrast_Callback(hObject, eventdata, handles)
    global postprocess;

    value = get(hObject, 'value');

    range = 0.5 + value;  
    imin = double(min(postprocess(:)));
    imax = double(max(postprocess(:)));
    postprocess = (postprocess - imin) / (imax - imin) * range;
    postprocess = min(1, max(0, postprocess));

    axes(handles.display);
    imshow(postprocess);

% --- Executes during object creation, after setting all properties.
function slider_contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
% --- Increases/decreases lightness of the image by a constant
function slider_lightness_Callback(hObject, eventdata, handles)
    global postprocess;

    value = get(hObject, 'value');

    postprocess = postprocess + (value - 0.5);
    postprocess = min(1, max(0, postprocess));

    axes(handles.display);
    imshow(postprocess);   
    

% --- Executes during object creation, after setting all properties.
function slider_lightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_lightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in b_gamma.
% --- Performs automatic gamma correction
function b_gamma_Callback(hObject, eventdata, handles)
    global postprocess;
    
    gray = rgb2gray(postprocess);
    grayscale = 0.25/mean(gray(:));
    postprocess = min(1, postprocess * grayscale);
    postprocess = srgb_gamma(postprocess);                                                         

    axes(handles.display);                
    imshow(postprocess);

% --- Executes on button press in b_he.
% --- Performs histrogram equalization on YCbCr color space
function b_he_Callback(hObject, eventdata, handles)
    global postprocess;

    ycbcr = rgb2ycbcr(postprocess);
    y = ycbcr(:,:,1);
    y = histeq(y);
    ycbcr(:,:,1) = y;

    postprocess = ycbcr2rgb(ycbcr);
    axes(handles.display);
    imshow(postprocess);


% --- Executes on button press in b_restore.
% --- Restores the image after corrections to sRGB conversion state
function b_restore_Callback(hObject, eventdata, handles)
    global srgb postprocess;

    set(handles.slider_saturation, 'value', 0.5);
    set(handles.slider_contrast, 'value', 0.5);
    set(handles.slider_lightness, 'value', 0.5);
    
    set(handles.rb_srgb, 'Value', 1);
    
    postprocess = srgb;
    axes(handles.display);
    imshow(postprocess);


% --- Executes on button press in b_lrotate.
function b_lrotate_Callback(hObject, eventdata, handles)
    global postprocess;

    postprocess = imrotate(postprocess, 90, 'bilinear');
    
    axes(handles.display);
    imshow(postprocess);


% --- Executes on button press in b_rrotate.
function b_rrotate_Callback(hObject, eventdata, handles)
    global postprocess;

    postprocess = imrotate(postprocess, 270, 'bilinear');        
    
    axes(handles.display);
    imshow(postprocess);


% --- Executes on button press in b_save.
% --- Saves image
function b_save_Callback(hObject, eventdata, handles)
    global postprocess;

    imshow(postprocess);
    imsave;
