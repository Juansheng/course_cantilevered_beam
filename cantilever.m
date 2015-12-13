clear all;

# beam params
# Alumnum beam with rect cross-section: 2cm x 2cm
global beam_len = 2;
global mass = 2700 * 0.02^2;
global EI = 70*10^9 * 0.02^4/12;

# FEM
global n_dof = 2; # number of DOFs: v, v'
global n_elem = 8; # number of elements

# elemental stiffness/mass matrix 
l = beam_len/n_elem;
stff_elem = [   12,    6*l,   -12,    6*l;
               6*l,  4*l^2,  -6*l,  2*l^2;
               -12,   -6*l,    12,   -6*l;
               6*l,  2*l^2,  -6*l,  4*l^2] * EI/l^3;               
mass_elem = [  156,   22*l,    54,  -13*l;
              22*l,  4*l^2,  13*l, -3*l^2;
                54,   13*l,   156,  -22*l;
             -13*l, -3*l^2, -22*l,  4*l^2] * mass*l/420;

# assembly
n_global = (n_elem + 1)*n_dof;
# declare global stiffness/mass matrix
g_stff = g_mass = zeros(n_global, n_global);
delta = 2*n_dof - 1;
for i = 1:n_elem
    st = 1 + (i-1)*n_dof;
    ed = st + delta;
    g_stff(st:ed, st:ed) += stff_elem;
    g_mass(st:ed, st:ed) += mass_elem;
end

# boundary condition: cantilevered beam
global n_trim = n_global - 2;
global e_stff = e_mass = zeros(n_trim, n_trim); # global for ODE solver
e_stff = g_stff(3:end, 3:end);
e_mass = g_mass(3:end, 3:end);

# eigen system
[modal, lamda] = eig(inv(e_mass)*e_stff); # NOTE: lamda unsorted!
printf("\n");
printf("Eigenvalues by FEM: \n\n");
lamda1 = lamda2 = zeros(n_trim, 1);
for i = 1:n_trim
    lamda1(i) = lamda(i,i);
end
[lamda_sorted, lamda_index] = sort(lamda1);
sqrt(lamda_sorted(1:5))

# exact solution
function res = lamda_func(n)
    global beam_len;
    global mass;
    global EI;
    res = ((2*n-1)/2*pi/beam_len)^2*sqrt(EI/mass);
end
printf("\n");
printf("Exact solution: \n\n");
k1 = 1.875;
lamda2(1) = (k1/beam_len)^2*sqrt(EI/mass);
k2 = 4.694;
lamda2(2) = (k2/beam_len)^2*sqrt(EI/mass);
for i = 3:n_trim
    lamda2(i) = lamda_func(i);
end
# print freqs
lamda2(1:5)

# modal shape
n_order = 3; # check the 3rd modal shape for now
modal_picked = zeros(n_trim/2+1, 1);
modal_picked(1) = 0; # modal fixed on the first node
# pick v and drop v'
modal_picked(2:end) = modal(1:2:end, lamda_index(n_order));
# nodes
beam_loc = 0:beam_len/n_elem:beam_len;
# show modal shape
figure(1);
plot(beam_loc, modal_picked, '-*b');

# ODE solver
global f = 80;
global omg = 5; # rad/s
global e_force = zeros(n_trim, 1);
global i_mass = inv(e_mass);

function yt = func(y, t)
    global n_trim;
    global i_mass;
    global e_stff;
    global f;
    global omg;
    global e_force;

    # load applied on the free end
    e_force(end-1) = -f*sin(omg*t);

    # ydot evaluation
    yt = zeros(2*n_trim, 1);
    yt(1:n_trim) = y(n_trim+1:end);
    yt(n_trim+1:end) = i_mass*(e_force - e_stff*y(1:n_trim));
end

# solving the equation
t = 0:0.02:2*(2*pi/omg);
y = lsode("func", zeros(2*n_trim, 1), t);

# plot v at the free end
figure(2);
plot(t, y(:,n_trim-1), '-*b')
