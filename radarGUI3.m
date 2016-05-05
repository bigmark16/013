
function varargout = radarGUI3(varargin)
% RADARGUI3 MATLAB code for radarGUI3.fig
%      RADARGUI3, by itself, creates a new RADARGUI3 or raises the existing
%      singleton*.
%
%      H = RADARGUI3 returns the handle to a new RADARGUI3 or the handle to
%      the existing singleton*.
%
%      RADARGUI3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RADARGUI3.M with the given input arguments.
%
%      RADARGUI3('Property','Value',...) creates a new RADARGUI3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before radarGUI3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to radarGUI3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help radarGUI3

% Last Modified by GUIDE v2.5 04-May-2016 22:12:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @radarGUI3_OpeningFcn, ...
                   'gui_OutputFcn',  @radarGUI3_OutputFcn, ...
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


% --- Executes just before radarGUI3 is made visible.
function radarGUI3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to radarGUI3 (see VARARGIN)

% Choose default command line output for radarGUI3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes radarGUI3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = radarGUI3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [w] = hann_window(N)
% create a hann (cosine squared) window
w = .5 + .5*cos(2*pi*((0:N-1).'/(N-1) - .5));

% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% ---- setup constants and parameters ----
c = 299e6; % (m/s) speed of light
cpi = 0.50; % (s) coherent processing interval 
fc = 2400e6; % (Hz) Center frequency (connect VCO Vtune to +5)
maxSpeed = 30; % (m/s) maximum speed to display
overlapFactor = 8; % (unitless) number of overlapped pulse windows (1 for no overlap)
ovsDop = 4; % (unitless) oversample factor for Doppler axis
Fs = 44100
% ----- end constants and parameters -----
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1
myRecording = [];
while get(hObject,'Value')
        t = 1;
        recObj = audiorecorder(44100,24,2);
        disp('step recording')
        recordblocking(recObj, t);
        disp('step processing');

        % Store data in double-precision array.

        myRecording = [myRecording; getaudiodata(recObj)];

        % read the raw myRecording data
        fprintf('Loading step recording...\n');
        [Y] = myRecording;

        % derived parameters
        N = round(cpi*Fs); % # of samples per pulse

        % the input appears to be inverted
        x = -Y(:,1);
        clear Y;

        % grab an integer number of overlapped frames
        M = floor(numel(x) / N * overlapFactor) - (overlapFactor) + 1;

        % compute axes parameters for the plot
        % note: the Doppler data is oversampled by ovsDop
        delta_f = (0:ovsDop*N/2-1).' / (ovsDop*N) * Fs; % Doppler freq. (Hz)
        lambda = c / fc; % wavelength (m)
        speed = delta_f * lambda / 2; % Doppler -> speed
        time = (1:M) * cpi / overlapFactor; % collection time (sec)

        % limit the speed axis to a reasonable range
        speed = speed(speed <= maxSpeed);
        nSpeed = numel(speed);

        % compute a Doppler window
        dopWin = hann_window(N);

        % compute the Doppler vs. time plot
        fprintf('Processing...\n');
        dti = zeros(nSpeed, M);
        for mIdx = 1:M
            xInds = (1:N).' + (mIdx-1)*floor(N/overlapFactor); % compute indices
            tmp = double(x(xInds)) .* dopWin; % apply Doppler window
            tmp = tmp - mean(tmp); % remove DC component if it exists
            tmp = fft(tmp, ovsDop*N); % compute oversampled Doppler spectrum
            dti(:,mIdx) = 20*log10( abs(tmp(1:nSpeed)) ); % grab result in dB
        end
        clear x;
        dti = dti.';

        % make Doppler vs. time plot
        %figure;        
        imagesc(speed,time,dti);
        colormap(jet(256));
        caxis(max(dti(:)) + [-60 0]); % show 60 dB dynamic range
        colorbar;
        xlabel('Speed (m/sec)');
        ylabel('Time (sec)');
        fprintf('step');0
end