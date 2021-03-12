classdef PcEllipsoidConductivity
    properties
        center
        radii
        axes
        amplitude
        n
    end
    methods
        function obj = PcEllipsoidConductivity(cen,radi,ax,amp,n)
            obj.center = cen;
            obj.radii = radi;
            obj.axes = ax;
            obj.amplitude = amp;
            obj.n = n;
        end
    end
end