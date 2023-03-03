function plot_juliareach(vars,x1,x2,clr,varargin)
%PLOT_JULIAREACH
if nargin > 4
    if strcmp(string(varargin{1}),"last") 
        sz = length(vars.timeZ);
        centers = vars.centers;
        gens = vars.gens;
        for i = sz:sz
            tg = reshape(gens(i,:,:),[size(gens,2),size(gens,3)]);
            hold on;
            Zz = Zono(centers(i,:)', tg);
            Ss = Zz.toStar;
            Star.plotBoxes_2D_noFill(Ss,x1,x2,clr);
        end
    else
        error('Wrong inputs')
    end
else
    sz = length(vars.timeZ);
    centers = vars.centers;
    gens = vars.gens;
    for i = 1:sz
        tg = reshape(gens(i,:,:),[size(gens,2),size(gens,3)]);
        hold on;
        Zz = Zono(centers(i,:)', tg);
        Ss = Zz.toStar;
        Star.plotBoxes_2D_noFill(Ss,x1,x2,clr);
    end
end

