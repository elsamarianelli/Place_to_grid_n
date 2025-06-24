%% 11/10/2016
% DORDEKS implementation of fast iterative shrinkage-thresholding algorithm
%% em COMMMENTS
% % RUNFISTA - Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) for NNPCA-like decomposition
%
% OVERVIEW:
% This function finds a set of spatial activity maps (e.g., grid-like patterns) 
% that are smooth, non-negative, and orthogonal — similar to Non-Negative PCA (NNPCA).
% It's based on the FISTA optimization method, and applies Difference of Gaussians (DoG) 
% to create center-surround receptive fields.
%
% STEP-BY-STEP ALGORITHM:
%
% For each desired component (e.g., grid cell):
%
% 1. Initialize a random non-negative spatial map (x).
% 2. Normalize x to have unit energy (||x|| = 1).
% 3. Set y = x (used for momentum step).
%
% 4. Repeat for a number of iterations:
%    a. Compute Difference of Gaussians (DoG) on y to highlight local structure.
%    b. Use that to calculate a gradient step.
%    c. Orthogonalize the gradient with respect to previously found components.
%    d. Apply the gradient step to y to get new candidate map q.
%    e. Enforce non-negativity on q (set negative values to zero).
%    f. Normalize q to unit norm.
%    g. Use FISTA momentum to update y based on current and previous x.
%    h. Optionally stop early if energy stabilizes.
%
% 5. Save the final map x as the new component.
% 6. Repeat for all components (e.g., grid cells).
%
% INPUTS:
%   dims         : [height, width] of spatial map
%   sigma_vector : [sigma_center, sigma_surround] for DoG filter
%   nCells       : number of components to extract
%   boundaries   : logical (true = periodic, false = zero-padding)
%   iterations   : max number of FISTA iterations per component
%               'replicate' or false 	Edges are padded with zeros – 
%                non-periodic environment (e.g., square or bounded space).
%               'circular' or true Edges wrap around – periodic space 
%               (e.g., a torus or donut).
% OUTPUT:
%   a            : [height x width x nCells] grid-like activity maps
%   Energy_array : energy per iteration for each component
            % higher energy means the current map (or component) contains strong, 
            % structured spatial features, like blobs or grid patterns.
            % FISTA tries to maximize this energy, while keeping the result smooth, orthogonal, and non-negative.
function [a, Energy_array ] = RunFISTA(dims,sigma_vector,nCells,boundaries,iterations) %#ok
%% inputs:
% dims : size of surrounding box
% sigma_vector : standard deviation of Gaussian blurs
% nCells : number of grid cell
% boundaries : do we have periodic boundary conditions?
% iterations: maximal number of FISTA iterations
%% Outputs
% x - The infered neuronal activity  [dims, nModules]
% MSE_array - mse error (residual): ||x-data||^2, at each itetation of FISTA

%  internal params
plotProgress=1;
tol=1e-20; %tolerance for stopping convergence
projectGrad=0;%should we project gradient on the tanget manifold to the sphere? 
noise=0;

D=length(dims);
Energy_array=zeros(iterations,nCells);
a=zeros([dims,nCells]);


% precompute some stuff
L=4; %Liphshitz constant
eta=2/L;% learning rate
win_size=max(sigma_vector)*5; %Gaussian window width
gauss_diff_conv=@(z) imblur(z,sigma_vector(1)*ones(1,D),win_size,D,boundaries) ...
            - imblur(z,sigma_vector(2)*ones(1,D),win_size,D,boundaries) ;

if plotProgress==1
    fh=figure(192);
    set(fh,'units','normalized','outerposition',[0.05 0.1 0.9 0.8]);
end
  
tic
%main loop 
for cc=1:nCells
    
    % initialize
    if D>1
        x=randn(dims);
    else
        x=randn(dims,1);
    end
    x(x<=0)=0;

    x=x./sqrt(sum(x(:).^2));
    t_next=1;
    y=x;

    % FISTA loop

    for kk=1:iterations
        %initialize             
            t=t_next;
            x_prev=x;

        %gradient update
        z=y;
        for nn=1:(cc-1)
            temp=a(:,:,nn).*z;
            z=z-a(:,:,nn)*sum(temp(:));
        end
        
        r=gauss_diff_conv(z);
        grad_step=(eta/cc)*gauss_diff_conv(r);

        for nn=1:(cc-1)
            temp=a(:,:,nn).*grad_step;
            grad_step=grad_step-a(:,:,nn)*sum(temp(:));
        end

        if projectGrad==1 %project gradient on the tanget manifold to the sphere
                grad_step=grad_step-bsxfun(@times,y,grad_step.*y); %project gradient on the tanget manifold to the sphere
        end
            q=y+grad_step;

            % non_negative
            q(q<0)=0;   
            % normaliation
            x=q./sqrt(sum(q(:).^2));


        % Some more Fista stuff
            t_next=(1+sqrt(1+4*t^2))/2;
            y=x+((t-1)/t_next)*(x-x_prev)+noise*randn(size(x));


        % Adaptive restart from Donoghue2012    
    %     do_restart=(mean2((q-y).*(x-x_prev)))>0;
    %     if do_restart
    %         t_next=1;
    %         y=x;        
    %     end
        current_energy=sum(r(:).^2);
        Energy_array(kk,cc)=current_energy;
        if kk>1
            relative_error=abs(1-Energy_array(kk,cc)/Energy_array(kk-1,cc));
            if relative_error<tol
                Energy_array(kk+1:end,cc)=nan;
                break
            end
        end
        if plotProgress==1;
            if ~mod(kk,50)
                rows=round(sqrt(nCells));
                col=ceil(nCells/rows);
                subplot(rows,col,cc)
                imagesc(squeeze(x(:,:)));            
    %             imagesc(squeeze(x(:,:,1+mod(kk,size(x,3)))));
                title(['Cell#' num2str(cc) ', Iter#' num2str(kk) ])
                pause(1e-4)
            end
        end
        if ~mod(kk,100)
            disp([ 'Cell #' num2str(cc) ', Iter #' num2str(kk)])
        end
    end
    if nCells>1
        a(:,:,cc)=x;
%         A=reshape(a(:,:,1:cc),[prod(dims),cc]);
    else
        a=x;
    end
end

toc