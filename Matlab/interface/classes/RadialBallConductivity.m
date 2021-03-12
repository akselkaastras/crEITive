classdef RadialBallConductivity
    properties
        radius
        amplitude
        innerextent
    end
    methods
        function obj = RadialBallConductivity(rad,amp,inext)
            obj.radius = rad;
            obj.amplitude = amp;
            obj.innerextent = inext;
        end
    end
end