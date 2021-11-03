function progress(count,N,step)
%
% progress(count,N,step)
%
% Placed within a loop, a call to the PROGRESS function displays a progress
% message and some time estimate. This estimate is only exact if each run
% of the loop takes roughly the same amount of time.
%
% Prior to the loop, call the function without arguments in order to
% initialize the function:
%
% progress();
%
% At the beginning of the body of the loop, call 
%
% progress(count,N,step);
%
% The variable COUNT specifies the number of the current run, 
% and N denotes the total number of runs of this loop. STEP is
% optional and determines that a progress message is issued
% at each STEPth run of the loop. Standard value for STEP is 20.
%
% Note: This function does not affect the use of TIC and TOC. 
%       Use TIC and TOC as usual in your code.

% initialize and retrieve persistent variable
persistent tmarker;

if nargin == 0 % we just need to initialize the timer and the starting time
    tmarker = tic;
    return
end

% some error checking
if nargin < 2
    error('Input arguments not defined.')
end

% assign standard values if variable is not defined
if nargin < 3
    step = 20;
end

% if we are not at the correct iteration step, we do not need 
% to report anything
if rem(count,step) ~= 0
    return
end

% if this is the first iteration, we still need to collect data in order to
% estimate the remaining amount of time
if count == 1
    fprintf(' Iteration: %i\t Collecting data in order to estimate time to completion. Computation startet at: %s\n',count,datestr(now))
    return
end

% If we reach this line here, we do our usual reporting.
% (count-1) instead of count is used because we expect the function to be
% called at the beginning of a loop and thus can only estimate the
% remaining time to completion for the following iteration.
sec = (N-(count-1))*(toc(tmarker)/(count-1));
fprintf(' Iteration: %i\t Estimated time to completion: %.3f minutes (= %f hours);\t Estimated end: %s\n',...
        count,sec/60,sec/3600,datestr(now+(sec/86400)))
