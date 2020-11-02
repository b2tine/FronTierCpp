clc
clear


%Fluid Data

timefile = dir('time*.txt');
meshfile = dir('mesh*.txt');
velfile = dir('velmat*.txt');
vortfile = dir('vortmat*.txt');

tdata = load(timefile.name);
mdata = load(meshfile.name);
veldata = load(velfile.name);
vortdata = load(vortfile.name);

string = split(velfile.name,'.');
string = split(string(1),'-');
size = cellfun(@str2num,string(2:4));

%TODO: add input option for specifying the non-periodic boundary dir
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

save('FluidData.mat','timeVec','mdata','velMat','vortMat');


%Interface Data

intfcposfile = dir('posintfc*.txt');
intfcvelfile = dir('velintfc*.txt');
intfcvortfile = dir('vortintfc*.txt');

iposdata = load(intfcposfile.name);
iveldata = load(intfcvelfile.name);
ivortdata = load(intfcvortfile.name);

string = split(intfcposfile.name,'.');
string = split(string(1),'-');
size = cellfun(@str2num,string(2:3));

stride = size(1);
tsteps = size(2);

%TODO: this won't work when stride is not a constant
iposMat = zeros(stride,2,tsteps);
ivelMat = zeros(stride,2,tsteps);
ivortMat = zeros(stride,tsteps);

for t = 1:tsteps
    iposMat(:,:,t) = iposdata((t-1)*stride+1:t*stride,:);
    ivelMat(:,:,t) = iveldata((t-1)*stride+1:t*stride,:);
    ivortMat(:,t) = ivortdata((t-1)*stride+1:t*stride);
end

save('IntfcData.mat','iposMat','ivelMat','ivortMat');

