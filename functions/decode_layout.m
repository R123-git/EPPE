function S=decode_layout(x,p)
assert(numel(x)==numel(p.available) && all(isfinite(x)));
[~,order]=sort(x(:),'ascend');
S=sort(p.available(order(1:p.k)));
end
