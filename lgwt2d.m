function [xx,yy,ww] = lgwt2d(x,w)
    [X,Y]=meshgrid(x,x);
    [W1,W2] = meshgrid(w,w);
    ww = reshape(W1.*W2,[numel(W1),1]);
    xx = reshape(X,[numel(X),1]);     yy = reshape(Y,[numel(Y),1]);
end
