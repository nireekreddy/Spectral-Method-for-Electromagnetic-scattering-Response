classdef CylinderGeometry
  properties
    % simulation properties
    orders@double;

    % background properties
    ep@double scalar complex;
    mu@double scalar complex;

    % field properties
    k@double scalar complex;
    beta@double scalar complex;
    phi@double scalar complex;

    % cylinder properties
    a@double scalar;
    epi@double scalar complex;
    mui@double scalar complex;

    % cylinder coordinates
    x@double scalar; 
    y@double scalar;

    % cylinder coefficients
    AmE; BmE; CmE;
    AmH; BmH; CmH;
  end

  methods
    % index
    function val = n(o)
      val = sqrt(o.ep*o.mu);
    end

    function val = ni(o)
      val = sqrt(o.epi*o.mui);
    end

    % in-plane propagation constant
    function val = kperp(o)
      val = sqrt([o.k].^2.*[o.ep].*[o.mu] - [o.beta].^2);
    end

    function val = kperpi(o)
      val = sqrt([o.k].^2.*[o.epi].*[o.mui] - [o.beta].^2);
    end

    % cylinder coordinates polar
    function rc = r(o)
      [~, rc] = cart2pol([o.x], [o.y]);
    end

    function thc = th(o)
      [thc, ~] = cart2pol([o.x], [o.y]);
    end

    % initialize empty coefficients
    function o = initcoeff(o)
      if isempty(o.AmE), o.AmE = zeros(size(o.orders)); end;
      if isempty(o.BmE), o.BmE = zeros(size(o.orders)); end;
      if isempty(o.CmE), o.CmE = zeros(size(o.orders)); end;
      if isempty(o.AmH), o.AmH = zeros(size(o.orders)); end;
      if isempty(o.BmH), o.BmH = zeros(size(o.orders)); end;
      if isempty(o.CmH), o.CmH = zeros(size(o.orders)); end;
    end
  end
end
