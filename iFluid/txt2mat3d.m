clc
clear

timefile = dir('time*.txt');
meshfile = dir('mesh*.txt');
velfile = dir('vel*.txt');
vortfile = dir('vort*.txt');

tdata = load(timefile.name);
mdata = load(meshfile.name);
veldata = load(velfile.name);
vortdata = load(vortfile.name);

string = split(velfile.name,'.');
string = split(string(1),'-');
size = cellfun(@str2num,string(2:5));

%TODO: add input option for specifying the non-periodic boundary dir
M = size(1) - 3;
N = size(2);
Q = size(3) - 3;
stride = M*N*Q;
tsteps = size(4);

timeVec = zeros(tsteps,1);
velMat = zeros(stride,3,tsteps);
vortMat = zeros(stride,3,tsteps);

for t = 1:tsteps
    timeVec(t) = tdata(t,1);
    velMat(:,:,t) = veldata((t-1)*stride+1:t*stride,:);
    vortMat(:,:,t) = vortdata((t-1)*stride+1:t*stride,:);
end

t_star = timeVec;
x_star = mdata(:,1);
y_star = mdata(:,2);
z_star = mdata(:,3);
u_star = zeros(stride,tsteps);
v_star = zeros(stride,tsteps);
w_star = zeros(stride,tsteps);
vortx_star = zeros(stride,tsteps);
vorty_star = zeros(stride,tsteps);
vortz_star = zeros(stride,tsteps);

for t = 1:tsteps
    u_star(:, t) = velMat(:,1,t);
    v_star(:, t) = velMat(:,2,t);
    w_star(:, t) = velMat(:,3,t);
    vortx_star(:, t) = vortMat(:,1,t);
    vorty_star(:, t) = vortMat(:,2,t);
    vortz_star(:, t) = vortMat(:,3,t);
end
 

save('FT_Cyl3d.mat','timeVec','mdata','velMat','vortMat');
save('FT_Cyl3d_sepcomp.mat', 't_star', 'x_star', 'y_star', 'z_star', ...
'u_star', 'v_star', 'w_star', 'vortx_star', 'vorty_star', 'vortz_star');
