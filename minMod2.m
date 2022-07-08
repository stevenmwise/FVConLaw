function m = minMod2(a,b)

m = (sign(a)==sign(b)).*sign(a).*min(abs([a b]),[],2);

end