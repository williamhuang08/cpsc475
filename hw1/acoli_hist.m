function [histm, foodfun] = ...
   acoli_hist(initpos, foodfn, rttp, ttrp, ntrials, nsteps, nbins, visualize)
%ACOLI_HIST Computes histograms for A. Coli 
%   [HISTM, FOOD] =
%     acoli_hist(INITPOS, FOODFN, RTTP, TTRP, NTRIALS, NSTEPS, NBINS, VIZ)
%
%   Computes the number of times a bacterium enters a region. The regions
%   are defined as rectangular bins which cover the space in a regular
%   grid. The experiment consists of several trials (ntrials) during each
%   of which the bacterium is allowed to move for several steps.
%   
%   Behavior: The bacterium either performs a run (advancing two units in
%   the current direction) or a tumble (rotating the current direction by a
%   random angle between 20 and 90 degrees and advancing one unit).
%   Switching between the two modes is controlled by the two probabilities
%   RTTP and TTRP. In addition, if the bacterium senses that it has moved
%   in the direction of increasing food concentration, then the random
%   angle during a tumble is between 10 and 45 degrees, i.e. the bacterium
%   does not turn much if it goes toward food.
%   
%   -- Input Arguments --
%
%   INITPOS is the initial position of the bacterium as a 2 element vector
%   [X Y]. The allowable range for both coordinates is [-2, 2]. To start
%   the bacterium in the center of the food function, use [0, 0].
%
%   FOODFN selects the food function. This can be either 1 or 2.
%
%   RTTP is the probability of transitioning from the run state to the
%   tumble state (a number between 0 and 1).
%
%   TTRP is the probability of transitioning from a tumble ot a run (a number 
%   between 0 and 1).
%
%   NTRIALS is the number of separate trials to run (starting the bacterium
%   from INITPOS each trial).
%
%   NSTEPS is the number of times to let the bacterium update its
%   behavior during each trial (i.e., the number of times to perform either
%   a run or a tumble).
%
%   NBINS are the number of bins to use in each direction for the
%   histogram. (The histogram covers the same area, so using more bins
%   increases the resolution of the histogram.)
%
%   VIZ is a flag controlling whether or not to open a figure window
%   displaying movements of the bacterium during a trial. This is optional
%   and defaults to false. Remember that you can use CTRL+C at any time to
%   abort a long running MATLAB operation, in case NTRIALS and/or
%   RUNSPERTRIAL are large. This display slows down the simulation
%   considerably, so don't use it if you only need the histogram.
%   
%   -- Output Variables --
%   
%   HISTM is the matrix of bins containing the number of times a bacterium
%   has passed through the corresponding region, i.e., the histogram
%   matrix. Each element contains the total number of times the bacterium
%   has passed through that bin over all trials. The matrix covers a 4x4
%   square (these are the same arbitrary units as passed to INITPOS),
%   centered at the origin.
%
%   FOOD is the matrix containing the concentration of food, i.e., the food
%   function. It covers the same area as the histogram matrix.

if ~exist('visualize', 'var')
  visualize = false;
end

recordhist = 1;
sp = 0.1;
xmin = -2; xmax = 2;
ymin = -2; ymax = 2;
[xs, ys] = meshgrid(xmin:sp:xmax, ymin:sp:ymax);

xbins = nbins;
ybins = nbins;

xlocs = linspace(xmin, xmax, xbins);
ylocs = linspace(ymin, ymax, ybins);

xstp = (xmax - xmin)/length(xlocs);
ystp = (ymax - ymin)/length(ylocs);

histm = zeros(xbins, ybins);

if foodfn == 1
  food.func = @foodfunc1;
elseif foodfn == 2
  food.func = @foodfunc2;
else
  error('foodfn must be either 1 or 2.') 
end
    
food.pts = [1 1; -0.5 -0.5];
[xx, yy] = meshgrid(linspace(xmin, xmax, xbins), linspace(ymin, ymax, ybins));
foodfun = food.func(xx, yy);

dt = sp;

for runi = 1:ntrials
  acoli.pos = initpos;
  acoli.oldpos = acoli.pos;
  acoli.vec = [rand - 0.5, rand - 0.5];
  acoli.fgrad = 0;
  acoli.mode = 0;
  acoli.hist = acoli.pos;

  for i = 1:nsteps
    acoli = update_bacterium(acoli, food, dt); 
    px = acoli.pos(1);
    py = acoli.pos(2);
    if px >= xmin && px <= xmax && py >= ymin && py <= ymax
      ptj = floor((px - xmin)/xstp) + 1;
      pti = floor((py - ymin)/ystp) + 1;
      histm(pti, ptj) = histm(pti, ptj) + 1;
    end

    if visualize
      contour(xs, ys, food.func(xs, ys));
      hold on;
      plot(food.pts(1,2), food.pts(1,2), '*');
      display_trajectory(acoli);
      colormap(0.7*[1 1 1]);
      hold off;
      drawnow;
    end
  end
end

%%
  function c = foodfunc1(x,y)
    s = 0.5;
    c = exp(-s*((x - food.pts(1,1)).^2 + (y - food.pts(1,2)).^2));
    
  end
%%
  function c = foodfunc2(x,y)
    s = 5;
    c = 2*exp(-s*((x - food.pts(1,1)).^2 + (y - food.pts(1,2)).^2)) + ...
        exp(-s*((x - food.pts(2,1)).^2 + (y - food.pts(2,2)).^2));
  end

%%
  function nv=rotatevec(v,a)
    rotM = [cos(a) sin(a); -sin(a) cos(a)];
    nv = (rotM*v(:))';
  end

%%
  function newacoli = update_bacterium(acoli, food, dt)
    fdist = food.func; % Food distribution function
    
    % Compute concentration gradient
    fgrad = fdist(acoli.pos(1), acoli.pos(2)) - ...
            fdist(acoli.oldpos(1), acoli.oldpos(2));
    
    newvec = acoli.vec;
    newacoli = acoli;
    
    rotupgrad = 1; % Scale factor for tumble rotation
    
    tdiff = food.func(acoli.pos(1), acoli.pos(2)) - ...
            food.func(acoli.oldpos(1), acoli.oldpos(2));
    
    if tdiff > 0
      rotupgrad =  2;
    end
    
    if acoli.mode == 0 % Run
      if rttp > rand
        newacoli.mode = 1; % Change to tumble mode
      end
      
      vec = acoli.vec/norm(acoli.vec);
      newvec = 2*vec; %Two units, same direction
      
    elseif acoli.mode == 1 % Tumble
      if ttrp > rand
        newacoli.mode = 0; % Change to run mode
      end
      
      vec = acoli.vec/norm(acoli.vec); % old direction
      rnum = rand;
      angle = ((1 - rnum)*20 + rnum*90)/rotupgrad*pi/180;
      newvec = rotatevec(vec, angle); % One unit, new direction
      
    end
    
    newpos = acoli.pos + dt*newvec; % Move
    
    newacoli.oldpos = acoli.pos;
    newacoli.pos = newpos;
    newacoli.vec = newvec;
    newacoli.fgrad = fgrad;
    
    if recordhist ~=0
      newacoli.hist = [acoli.hist; newpos];
    end
  end

%%
  function display_trajectory(acoli)
    % Draw trajectory
    plot(acoli.hist(:,1), acoli.hist(:,2), ':');
    
    % Draw initial location
    plot(acoli.hist(1,1), acoli.hist(1,2), 'ks', 'MarkerFaceColor', 'k');
    
    % Draw end location
    plot(acoli.pos(1), acoli.pos(2), 'ko', 'MarkerFaceColor', 'r');
    
    % Draw the direction of the acoli
    vec = acoli.vec;
    vec = vec/norm(vec)*.2;
    quiver(acoli.pos(1),acoli.pos(2),vec(1),vec(2));
    
    axis equal tight
  end
end
