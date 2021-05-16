% Class
classdef PcBallShell
    properties
        center
        radius
        amplitude
        width
        n
    end
    methods
        function obj = PcBallShell(cen,rad,amp,width,n)
            obj.center = cen;
            obj.radius = rad;
            obj.amplitude = amp;
            obj.width = width;
            obj.n = n;
        end
    end
end
