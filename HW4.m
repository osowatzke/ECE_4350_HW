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
fprintf("(1)\n\n");
T = table(modes',modes_fc','VariableNames',{'mode','fc (GHz)'});
disp(T);

%% problem 2

% waveguide dimensions in inches
a = 5.1;
b = 2.55;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% determine cutoff wavelength
Lco = 2*a;

% output cutoff wavelength
fprintf("2a) Cutoff Wavelength: %.4f m\n", Lco);

% permeability 
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% determine cutoff frequency
fco = 1/(2*a*sqrt(u*e));

% output cutoff frequency
fprintf("2a) Cutoff Frequency: %.4f GHz\n",fco*1e-9);

% frequency
f = 10e9;

% determine wavenumber
k = 2*pi*f*sqrt(u*e);

% determine propagation constant
B = k*sqrt(1-(fco/f)^2);

% output propation constant
fprintf("2b) Progation Constant: %.4f rad/m\n",B);

% determine wavelength
L = (1/sqrt(u*e))/f;

% determine characteristic impedance of dielectric
n = sqrt(u/e);

% Transverse-wave impedance
ZTE = n*(1-(L/(2*a))^2)^(-1/2);

% output transverse-wave impedance
fprintf("2b) Transverse-Wave Impedance: %.4f Ohms\n", ZTE);

% determine TEM propagation constant
B = k;

% determine TEM transverse wave impedance
ZTEM = n;

% output propogation constant
fprintf("2c) TEM Progation Constant: %.4f rad/m\n",B);

% output transverse wave impedance
fprintf("2c) TEM Transverse-Wave Impedance: %.4f Ohms\n",ZTEM);

%% problem 3








