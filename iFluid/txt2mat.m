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
size = cellfun(@str2num,string(2:4));

M = size(1);
N = size(2) - 3;
stride = M*N;
tsteps = size(3);

timeVec = zeros(tsteps,1);
velMat = zeros(stride,2,tsteps);
vortMat = zeros(stride,tsteps);

for t = 1:tsteps
    timeVec(t) = tdata(t,1);
    velMat(:,:,t) = veldata((t-1)*stride+1:t*stride,:);
    vortMat(:,t) = vortdata((t-1)*stride+1:t*stride);
end

save('FT_cyl2d.mat','timeVec','mdata','velMat','vortMat');
