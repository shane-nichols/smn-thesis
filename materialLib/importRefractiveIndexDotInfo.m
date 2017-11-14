% for isotropic stuff
f = fopen('Ag_raw.txt');

g = textscan(f,'%s','delimiter','\n');
g = squeeze(g{1}).';
fclose(f);
idx = 0;
for i = g
    idx = idx+1;
    if isempty(find(strcmpi(i,'wl	k'), 1)) == 0
        breakAt = idx;
    else
    end
end

for i=1:(breakAt-10)
    n(i,:) = str2num(g{i+1});
    k(i,:) = str2num(g{i+1+breakAt});
end

Ag=[n,k(:,2)];
save Ag.txt Ag -ascii