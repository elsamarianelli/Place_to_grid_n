function varargout = tintplot( varargin )
% See https://github.com/UCL/mTint/tree/master/LoadData for info

    if nargin == 0, varargout{1} = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    if ~isvector(in.im) % dealing with a ratemap
        % Get a colormap of the desired size (see subfunction)
        [cm,NUM_COLORS] = tintColorMap(in.NUM_COLORS);

        % Apply the colormap to the image
        im = in.im;
        bad = isnan(im);
        if isempty(in.CLIM)
            m = max(im(~bad(:)));
        else
            m = opts.CLIM;
            im(im>m) = m;
        end
        im = round(im/m*(NUM_COLORS+0.999) + 0.5);
        img = ind2rgb(im,cm);

        % Set non-visited pixels to be white
        bad = repmat(bad,[1 1 3]);
        img(bad) = 1; 

        % Display the image
        if ~isempty(in.AXIS_LIMS)
            %here the user specifies the cm value of the centres of the most
            %extreme bins in each dimensions
            image(in.AXIS_LIMS(1,:),in.AXIS_LIMS(2,:),img); 
            xlabel('cm'); ylabel('cm');
        else
            % here we don't bother with true units, and just use "bins"
            image(img);
            xlabel('bins'); ylabel('bins');
        end

        axis image ij
        
    else% dealing with a polar plot
        rho = in.im;
        nBins = length(rho);
        theta = 0:2*pi/nBins:2*pi-(2*pi)/nBins;
        polar(theta, rho')
        delete(findall(ancestor(gcf,'figure'),'HandleVisibility','off',...
            'type','line','-or','type','text'));
        ylim = get(gca,'ylim');
        xlim = get(gca,'xlim');
        xHair = min(abs([xlim,ylim])) * [-1 1];
        hold on
        vert_line = line([0,0],xHair);
        horz_line = line(xHair,[0,0]);
        set(vert_line,'Color','k')
        set(horz_line,'Color','k')
        tick_len = max(xHair) * 0.02;
        ticks = linspace(xHair(1), xHair(2), 21);
        vert_ticks = line([zeros(1,21)-tick_len;zeros(1,21)+tick_len], [ticks; ticks]);
        set(vert_ticks,'Color','k')
        horz_ticks = line([ticks; ticks],[zeros(1,21)-tick_len;zeros(1,21)+tick_len]);
        set(horz_ticks,'Color','k')
        maxStr = num2str(max(rho));
        str = text(max(xHair) * 0.05, ticks(end-1), [maxStr(1:strfind(maxStr,'.')+2),' Hz']);
        set(str, 'Fontsize', 8);
    end
end


function [cbjet,NUM_COLORS]=tintColorMap(NUM_COLORS)
%only 5 at the moment

switch NUM_COLORS
    case 5
         cbjet=[0,  0, 0.7765; ...
        0,  0.6353, 1; ...
        0.2196,  .9216, 0.1255; ...
        .972,  .86706, 0; ...
        1, .1255, 0]; %NB. This should be identical to Tint's 5-color setting
    otherwise
         cbjet=jet(NUM_COLORS);
         %{
            [0,  0, 0.7765; ...
                0,  0.6353, 1; ...
                0.2196,  .9216, 0.1255; ...
                .972,  .86706, 0; ...
                1, .1255, 0]; %NB. This should be identical to Tint's 5-color setting
          NUM_COLORS = 5;
         %}
end
    %If you want to do more than 5...
    %do an example plot in Tint, use printscreen to copy it.
    %paste in paint and use color picker to select color.
    %go to color>>edit colors>>define custom colors and look at R, G and B
    %divde these numbers by 255 and put them here in a matrix.
end


function in = DealWithInputs(varargin)
    defaults.im = []; 
    defaults.NUM_COLORS = 5; 
    defaults.AXIS_LIMS = []; 
    defaults.CLIM = [];
    
    VERSION = 1;

    if nargin == 0
        in = defaults;
        return;
    end

    if nargin > 1 && VERSION ~= varargin{2}
        error(['tintplot called with version number %g, but code is version %g.\n' ...
            'Check the GitHub wiki to see what has changed, or take a risk and call the function without a version number.'],varargin{2},VERSION);
    end

    in = setstructfields(defaults,varargin{1});
end
