% Class
classdef PcBallConductivity
    properties
        center
        radius
        amplitude
        n
    end
    methods
        function obj = PcBallConductivity(cen,rad,amp,n)
            obj.center = cen;
            obj.radius = rad;
            obj.amplitude = amp;
            obj.n = n;
        end
    end
end

