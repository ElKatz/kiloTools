function [xcoords, ycoords, kcoords] = probeGeometry2coords(probeGeometry, nCh)
%   [xcoords, ycoords, kcoords] = probeGeometry2coords(probeGeometry, nCh)
%
% for a given probe geometry (string), function returns the x y and k
% coordinates for kiloSorting.
% INPUT:
%   probeGeometry - e.g. 'linear50' or 'dualLinear'. see optoins in code.
%   nCh - number of channles on the probe

    
switch probeGeometry
    case 'linear50'
        % linear geometry:
        ySpacing = 50;
        xcoords = zeros(nCh,1);
        ycoords = (((nCh-1) * ySpacing):-ySpacing:0)';
        kcoords = ones(nCh,1);
    case 'linear100'
        % linear geometry:
        ySpacing = 100;
        xcoords = zeros(nCh,1);
        ycoords = (((nCh-1) * ySpacing):-ySpacing:0)';
        kcoords = ones(nCh,1);
    case 'linear125'
        % linear geometry:
        ySpacing = 125;
        xcoords = zeros(nCh,1);
        ycoords = (((nCh-1) * ySpacing):-ySpacing:0)';
        kcoords = ones(nCh,1);
    case 'linear200'
        % linear geometry:
        ySpacing = 200;
        xcoords = zeros(nCh,1);
        ycoords = (((nCh-1) * ySpacing):-ySpacing:0)';
        kcoords = ones(nCh,1);
    case 'stereo'
        % sterotrode geometry
        xSpacing = 50; 
        ySpacing = 50;
        xcoords = xSpacing * repmat([0 1]', nCh/2,1);
        ycoords = ySpacing * kron(((nCh/2):-1:1)-1, [1 1])'; % using kron for stereo geometry
        kcoords = ones(nCh,1);
    case 'dualLinear'
        % two independent probes:
        nCh1 = nCh/2;
        nCh2 = nCh/2;
        xSpacing = 200; 
        ySpacing = 50;
        xcoords = xSpacing * [repmat(0, nCh1,1);            repmat(1, nCh2,1)];
        ycoords = [(((nCh1-1) * ySpacing):-ySpacing:0)';  (((nCh2-1) * ySpacing):-ySpacing:0)'];
        kcoords = ones(nCh,1);
    case 'dualLinearMinus1'
        % two independent probes (with the last channel accidentally
        % deactivated...)
        nCh1 = ceil(nCh/2);
        nCh2 = floor(nCh/2);
        xSpacing = 200; 
        ySpacing = 50;
        xcoords = xSpacing * [repmat(0, nCh1,1);            repmat(1, nCh2,1)];
        ycoords = [(((nCh1-1) * ySpacing):-ySpacing:0)';  (((nCh2-1) * ySpacing):-ySpacing:0)'];
        kcoords = ones(nCh,1);
    case 'vprobe50 plus injectrode' % GC added on 12/01/19 for vprobe inactivation experiment in SC
        % the channel for the vprobe is from 1 to (nCh-1) 
        % the 'nCh' th channel is the injectrode
        nCh1 = nCh - 1;
        nCh2 = 1;
        xSpacing = 200; 
        ySpacing = 50;
        xcoords = xSpacing * [repmat(0, nCh1,1);            repmat(1, nCh2,1)];
        ycoords = [(((nCh1-1) * ySpacing):-ySpacing:0)';  (((nCh2-1) * ySpacing):-ySpacing:0)'];
        kcoords = ones(nCh,1);
    otherwise
        error('Supplied probe geometry does not correspond to any geometry in this micro-universe!')
end