%% problem 1

% waveguide dimensions (in)
a = 0.9;
b = 0.4;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% number of modes desired
num_modes = 10;

% create row vector of all m and n values
m = zeros(1,(num_modes+1)^2);
n = zeros(1,(num_modes+1)^2);

% loop initializes row vectors
count = 1;
for i = 0:num_modes
    for j = 0:num_modes
        n(count) = i;
        m(count) = j;
        count = count + 1;
    end
end

% permeability
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% determine vector of cutoff frequencies
fc = 1./(2.*pi.*sqrt(u.*e)).*((m.*pi./a).^2+(n.*pi./b).^2).^(1/2);

% sort frequencies and indices
[fc,I] = sort(fc);
m = m(I);
n = n(I);

% loop through indices to find first 10 modes
modes_found = 0;
index = 2;
modes = strings(1,num_modes);
modes_fc = zeros(1,num_modes);
while (modes_found < 10)
    modes_found = modes_found + 1;
    modes(modes_found) = strcat("TE",int2str(m(index)),...
            int2str(n(index))); 
    modes_fc(modes_found) = fc(index)*1e-9;
    if m(index) ~= 0 && n(index) ~= 0
        modes_found = modes_found + 1;
        modes(modes_found) = strcat("TM",int2str(m(index)),...
            int2str(n(index)));
        modes_fc(modes_found) = fc(index)*1e-9;
    end
    index = index+1;
end

% output table of modes
T = table(modes',modes_fc','VariableNames',{'mode','fc (GHz)'});
disp(T);






