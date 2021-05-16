classdef PcEllipsoidShell
    properties
        center
        radii
        axes
        amplitude
        width
        n
    end
    methods
        function obj = PcEllipsoidConductivity(cen,radi,ax,amp,width,n)
            obj.center = cen;
            obj.radii = radi;
            obj.axes = ax;
            obj.amplitude = amp;
            obj.width = width;
            obj.n = n;
        end
    end
end