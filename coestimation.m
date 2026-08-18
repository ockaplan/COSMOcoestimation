% coestimation.m
% by Ollie Kaplan
% 7/22/26

%co-estimating magnetic field expansion coefficients alongside magnetometer
%mounting euler angles using Swarm L1b data via the method used by Olsen et
%al in 2000 for Orsted.

%this script REQUIRES:
%Parallelization Toolbox
%Aerospace Toolbox

clear; clc; close all;
rng(31770);

%% DATA INGEST

tic;

%get Swarm data
[Timestamp,Latitude,Longitude,Radius,F,B_VFM,B_NEC,q_CRF_NEC,q_NEC_CRF,eul_true,q_CRF_VFM] = ...
    getDataSet(60,"SW_OPER_MAGA_LR_1B_20260702T000000_20260702T235959_0702_MDR_MAG_LR.cdf");%,...
    %"SW_OPER_MAGA_LR_1B_20260703T000000_20260703T235959_0702_MDR_MAG_LR.cdf",...
    %"SW_OPER_MAGA_LR_1B_20260718T000000_20260718T235959_0702_MDR_MAG_LR.cdf");

%get Dst info by interpolation, data source: https://wdc.kugi.kyoto-u.ac.jp/dst_realtime/presentmonth/index.html

% 7/2 and 7/3
Dst_kyoto_23 = [7   6   7   6   6   5   2   2   3    3   5   7   7   4   0  -2  -2   -2   1   6   7   7   4   2   2 ...%07/02
    5   4   4   4   3   4   3   3    4   3   6   8  15  10   9   8    7  12   7   7   7   7  14  23  33]'; %07/03
time_kyoto_23 = datetime(2026,7,2,0,0,0) + hours(-0.5:48.5);

%7/18
Dst_kyoto_18 = [-6  -3  -1   2   2   5   7   5   3    4   6   7   7   4   0  -3  -3   -1   0   2   3   2   1   1   1   2]';
time_kyoto_18 = datetime(2026,7,18,0,0,0) + hours(-0.5:24.5);

Dst_kyoto = [Dst_kyoto_23;Dst_kyoto_18];
time_kyoto = [time_kyoto_23 time_kyoto_18];
Dst = interp1(time_kyoto,Dst_kyoto,Timestamp,"makima","extrap");

%load WMM2020 coefficients
coefMat = readmatrix("WMM2020COF/WMM.COF",'FileType','text');
coefMat = coefMat(1:90,:); %trim NaNs

%load WMM2025 coefficients
coefMat2026 = readmatrix("WMM2025.COF",'FileType','text');
coefMat2026 = coefMat2026(1:90,:); %trim NaNs

%start logging
logName = "logs/log" + string(datetime('now', 'Format', 'yyyyMMddHHmmss')) + ".txt";
logFile = fopen(logName,'w');
fprintf(logFile,"%s: Dataset loaded and trimmed.\n",strCurTime());

%% COESTIMATION

%overall estimation parameters
numIterations = 2;
maxDeg = 10;
numModelTerms = maxDeg*(maxDeg+2) + 10; %main field coefs + 1 int. Dst + 3 ext. + 3 ext Dst + 3 EAs
numMeas = height(Longitude);

%set up initial guess via WMM2020
dTrue = formTrueDataVec(F,B_VFM);
eul0 = eul_true + arcsec2rad(3600*rand(3,1) - 1800);
m = zeros(numModelTerms,numIterations+1);
initModel = getInitModelVec(coefMat,eul0,maxDeg);
m(:,1) = initModel.m;
fprintf(logFile,"%s: Formed initial model.\n",strCurTime());

%get what the model vec would be from WMM2025
m_wmm2026 = getInitModelVec(coefMat2026,eul_true,maxDeg);

%set up weight matrix
sigmaSc = 0.06; %nT
sigmaVec = 0.25; %nT
W = getWeightMatrix(Latitude,sigmaSc,sigmaVec);
fprintf(logFile,"%s: Formed weight matrix.\n",strCurTime());

%get terms for transforming coefficients -> VFM measurements
termMat = getTermMat(Longitude,Latitude,Radius,maxDeg,numModelTerms,Dst);
fprintf(logFile,"%s: Formed spherical harmonic expansion term matrix.\n",strCurTime());
attMatCRF = getMatCRF(q_NEC_CRF);
fprintf(logFile,"%s: Formed NEC->CRF transformation matrix.\n",strCurTime());

%iterate on the model
fprintf(logFile,"%s: Beginning iterative model refinement.\n",strCurTime());
for i = 1:numIterations
    modelVec = modelVecStruct(m(:,i),maxDeg);
    d_mi = getDataVec(termMat,attMatCRF,modelVec,numMeas);
    fprintf(logFile,"%s: Computed pre-update data vector for iteration %d. Beginning Jacobian.\n",strCurTime(),i);
    Gi = approxJacobian(modelVec.gVec,modelVec.eul,numMeas,termMat,attMatCRF,"FFD");
    fprintf(logFile,"%s: FFD Computed Jacobian for iteration %d.\n",strCurTime(),i);
    delta_mi = inv(Gi' * W * Gi) * (Gi' * W * (dTrue - d_mi));
    m(:,i+1) = m(:,i) + delta_mi;
    fprintf(logFile,"%s: Iteration %d complete.\n",strCurTime(),i);
end
totalElapsed = toc;
fprintf(logFile,"%s: Done.\n",strCurTime());

%get final predicted data vector
d_final = getDataVec_Vecinput(m(:,end),termMat,attMatCRF,numMeas);
resids = d_final-dTrue;

%plot evolution of EAs
figure(1)
subplot(3,1,1)
plot(m(end-2,:))
hold on
yline(eul_true(1))
subplot(3,1,2)
plot(m(end-1,:))
hold on
yline(eul_true(2))
subplot(3,1,3)
plot(m(end,:))
hold on
yline(eul_true(3))

%get euler angles between estimated mounting, initial estimate, and true
q_CRF_VFM_0 = angle2quat(eul0(1),eul0(2),eul0(3),'ZYX');
q_CRF_VFM_est = angle2quat(m(end-2,end),m(end-1,end),m(end,end),'ZYX');
[ang3_errInit, ang2_errInit, ang1_errInit] = quat2angle(quatmultiply(quatinv(q_CRF_VFM),q_CRF_VFM_0));
[ang3_errEst, ang2_errEst, ang1_errEst] = quat2angle(quatmultiply(quatinv(q_CRF_VFM),q_CRF_VFM_est));
eulErrorInit = rad2arcsec([ang3_errInit, ang2_errInit, ang1_errInit]);
eulErrorEst = rad2arcsec([ang3_errEst, ang2_errEst, ang1_errEst]);

%end logs
fprintf(logFile,"%s: Mounting quaternion residuals reduced from %.2f arcsec to %.2f arcsec.\n",strCurTime(),norm(eulErrorInit),norm(eulErrorEst));
fclose(logFile);

%save model
modelName = "model" + string(datetime('now', 'Format', 'yyyyMMddHHmmss')) + ".mat";
save("models/" + modelName,"m","maxDeg","numMeas","eulErrorInit","eulErrorEst","numIterations","totalElapsed","Timestamp")

%% PLOTS

F_resids = resids(1:numMeas);
X_resids = resids(numMeas+1:3:end);
Y_resids = resids(numMeas+2:3:end);
Z_resids = resids(numMeas+3:3:end);
f1 = figure();
subplot(4,1,1)
scatter(Latitude,F_resids,'Marker','.')
ylabel("F Residuals [nT]")
subplot(4,1,2)
scatter(Latitude,X_resids,'Marker','.')
ylabel("X Residuals [nT]")
subplot(4,1,3)
scatter(Latitude,Y_resids,'Marker','.')
ylabel("Y Residuals [nT]")
subplot(4,1,4)
scatter(Latitude,Z_resids,'Marker','.')
xlabel("Latitude [deg]")
ylabel("Z Residuals [nT]")
print("figures/Bcomps" + string(datetime('now', 'Format', 'yyyyMMddHHmmss')) + ".png","-dpng")

f2 = figure();
subplot(4,1,1)
scatter(Dst,F_resids,'Marker','.')
ylabel("F Residuals [nT]")
subplot(4,1,2)
scatter(Dst,X_resids,'Marker','.')
ylabel("X Residuals [nT]")
subplot(4,1,3)
scatter(Dst,Y_resids,'Marker','.')
ylabel("Y Residuals [nT]")
subplot(4,1,4)
scatter(Dst,Z_resids,'Marker','.')
xlabel("Dst [nT]")
ylabel("Z Residuals [nT]")
print("figures/DSTresids" + string(datetime('now', 'Format', 'yyyyMMddHHmmss')) + ".png","-dpng")

%% FUNCTION DEFS

function [Timestamp,Latitude,Longitude,Radius,F,B_VFM,B_NEC,q_CRF_NEC,q_NEC_CRF,eul_true,q_CRF_VFM] = getDataSet(meas2Skip,varargin)
    Timestamp = [];
    Latitude = [];
    Longitude = [];
    Radius = [];
    F = [];
    B_VFM = [];
    B_NEC = [];
    q_CRF_NEC = [];
    q_NEC_CRF = [];
    for i = 1:nargin-1
        %read cdf file
        fileStr = "data/" + varargin{i};
        %info = cdfinfo(fileStr);
        rawData = cdfread(fileStr,'DatetimeType','datetime');
        metaData1 = readstruct("data/SW_OPER_MAGA_LR_1B_20260702T000000_20260702T235959_0702.HDR", 'FileType', 'xml');
        
        %import data as MATLAB variables
        Timestamp_i = [rawData{:,1}]';
        %Timestamp.TimeZone = 'UTCLeapSeconds';
        %SyncStatus = [rawData{:,2}]';
        Latitude_i = [rawData{:,3}]';
        Longitude_i = [rawData{:,4}]';
        Radius_i = [rawData{:,5}]';
        F_i = [rawData{:,6}]';
        %dF_Sun = [rawData{:,7}]';
        %dF_AOCS = [rawData{:,8}]';
        %dF_other = [rawData{:,9}]';
        %F_error = [rawData{:,10}]';
        B_VFM_i = [rawData{:,11}]';
        B_NEC_i = [rawData{:,12}]';
        %dB_Sun = [rawData{:,13}]';
        %dB_AOCS = [rawData{:,14}]';
        %dB_other = [rawData{:,15}]';
        %B_error = [rawData{:,16}]';
        q_CRF_NEC_i = [rawData{:,17}]';
        q_CRF_NEC_i = [q_CRF_NEC_i(:,4) q_CRF_NEC_i(:,1:3)]; %change to scalar first
        q_NEC_CRF_i = quatinv(q_CRF_NEC_i);
        %Att_error = [rawData{:,18}]';
        %Flags_F = [rawData{:,19}]';
        Flags_B_i = [rawData{:,20}]';
        Flags_q_i = [rawData{:,21}]';
        Flags_Platform_i = [rawData{:,22}]';
        %ASM_Freq_Dev = [rawData{:,23}]';
        
        %get q_CRF_VFM from metadata
        q_VFM_CRF = [metaData1.Variable_Header.SPH.Magnetic_Information.q_STR_VFM.Q4,...
            metaData1.Variable_Header.SPH.Magnetic_Information.q_STR_VFM.Q1,...
            metaData1.Variable_Header.SPH.Magnetic_Information.q_STR_VFM.Q2,...
            metaData1.Variable_Header.SPH.Magnetic_Information.q_STR_VFM.Q3];
        q_CRF_VFM = quatinv(q_VFM_CRF);
        [ang3, ang2, ang1] = quat2angle(q_CRF_VFM,"ZYX");
        eul_true = [ang3, ang2, ang1]';
        
        %filter data with bad flags
        indB = Flags_B_i > 0;
        indq = Flags_q_i > 15;
        indPlat = Flags_Platform_i > 1;
        indLat = abs(Latitude_i) > 40; %filter out data from high latitudes
        ind = false(size(indB));
        for n=1:height(indB)
            ind(n) = ~(indB(n) || indq(n) || indPlat(n) || indLat(n));
        end
        Timestamp = [Timestamp;Timestamp_i(ind)];
        Latitude = [Latitude;Latitude_i(ind)];
        Longitude = [Longitude;Longitude_i(ind)];
        Radius = [Radius;Radius_i(ind)];
        F = [F;F_i(ind)];
        B_VFM = [B_VFM;B_VFM_i(ind,:)];
        B_NEC = [B_NEC;B_NEC_i(ind,:)];
        q_CRF_NEC = [q_CRF_NEC;q_CRF_NEC_i(ind,:)];
        q_NEC_CRF = [q_NEC_CRF;q_NEC_CRF_i(ind,:)];
    end
    %reduce number of measurements for computation time
    indTrim = 1:meas2Skip:height(Timestamp);
    Timestamp = Timestamp(indTrim);
    Latitude = Latitude(indTrim);
    Longitude = Longitude(indTrim);
    Radius = Radius(indTrim);
    F = F(indTrim);
    B_VFM = B_VFM(indTrim,:);
    B_NEC = B_NEC(indTrim,:);
    q_CRF_NEC = q_CRF_NEC(indTrim,:);
    q_NEC_CRF = q_NEC_CRF(indTrim,:);
end

function str = strCurTime()
    str = string(datetime("now"));
end

function dTrue = formTrueDataVec(F,B_VFM)
    %form measurements into true data vector
    B_VFM_t = B_VFM';
    dTrue = [F;B_VFM_t(:)];
end

function modelVec = getModelVec(gCoefs,g0Dst,hCoefs,qCoefs,qDst,sCoefs,sDst,eul)
    %convert model terms into modelVec struct
    modelVec.gnm = gCoefs; %(N^2+3N)/2 x 1
    modelVec.g0Dst = g0Dst; %1x1
    modelVec.hnm = hCoefs; %(N/2)(N+1) x 1
    modelVec.qnm = qCoefs; %2x1
    modelVec.qDst = qDst; %2x1
    modelVec.snm = sCoefs; %1x1
    modelVec.sDst = sDst; %1x1
    modelVec.eul = eul; %3x1
    modelVec.gVec = [gCoefs;hCoefs;qCoefs;sCoefs;g0Dst;qDst;sDst];
    modelVec.m = [modelVec.gVec;eul];
end

function modelVec = modelVecStruct(m,maxDeg)
    %convert m vector output by gauss estimator into modelVec struct
    numg = (maxDeg * maxDeg + 3 * maxDeg) * 0.5;
    numh = (maxDeg * maxDeg + maxDeg) * 0.5;
    gCoefs = m(1:numg,1);
    hCoefs = m(numg+1:numg+numh,1);
    qCoefs = m(end-9:end-8,1);
    sCoefs = m(end-7,1);
    g0Dst = m(end-6,1);
    qDst = m(end-5:end-4,1);
    sDst = m(end-3,1);
    eul = m(end-2:end,1);
    modelVec = getModelVec(gCoefs,g0Dst,hCoefs,qCoefs,qDst,sCoefs,sDst,eul);
end

function m0 = getInitModelVec(coefMat,eul0,maxDeg)
    %build initial model vec from WMM2020 and incorrect error angles
    endInd = 0.5*(maxDeg * maxDeg + 3 * maxDeg);
    maxIndWMM = min(endInd,0.5*(12 * 12 + 3*12));
    if maxDeg > 12
        coefMat = [coefMat; zeros(endInd-maxIndWMM,6)];
    end
    t = decyear(datetime(2026,7,2,0,0,0));
    g = coefMat(1:endInd,3) + (t - 2020) * coefMat(1:endInd,5);
    h = coefMat(1:endInd,4) + (t - 2020) * coefMat(1:endInd,6);
    if maxDeg > 12
        h([(coefMat(1:maxIndWMM,2) == 0);true;false(endInd-maxIndWMM-(maxDeg-12),1)]) = []; %remove terms for m=0
    else
        h((coefMat(1:maxIndWMM,2) == 0)) = [];
    end
    qCoefs = [0;0];
    sCoefs = 0;
    g0Dst = 0.23;
    qDst = [0.77;0];
    sDst = 0;
    m0 = getModelVec(g,g0Dst,h,qCoefs,qDst,sCoefs,sDst,eul0);
end

function block = getTBlock(Lon,Lat,r,Dst,maxDeg,numModelTerms)
    block = zeros(3,numModelTerms);
    a = 6371200; %mean magnetic radius [m]
    %populate g and h terms
    numg = (maxDeg * maxDeg + 3 * maxDeg) * 0.5;
    %numh = (maxDeg * maxDeg + maxDeg) * 0.5;
    gInd = 1;
    hInd = 1+numg;
    sinVec = sind(Lon * (0:maxDeg));
    cosVec = cosd(Lon * (0:maxDeg));
    secLat = secd(Lat);
    tanLat = tand(Lat);
    cosLat = cosd(Lat);
    sinLat = sind(Lat);
    for n = 1:maxDeg
        legVec = legendre(n,sind(Lat),'sch');
        legVecp1 = legendre(n+1,sind(Lat),'sch');
        for m = 0:n
            legDeriv = (n+1) * tanLat * legVec(m+1) - sqrt((n+1)^2 - m^2) * secLat * legVecp1(m+1);
            %g terms
            block(1,gInd) = -(a/r)^(n+2) * cosVec(m+1) * legDeriv;
            block(2,gInd) = 1/cosLat * (a/r)^(n+2) * m * sinVec(m+1) * legVec(m+1);
            block(3,gInd) = -(n+1) * (a/r)^(n+2) * legVec(m+1) * cosVec(m+1);
            gInd = gInd + 1;
            if m~=0
                %h terms
                block(1,hInd) = -(a/r)^(n+2) *  sinVec(m+1) * legDeriv;
                block(2,hInd) = 1/cosLat * (a/r)^(n+2) * -m * cosVec(m+1) * legVec(m+1);
                block(3,hInd) = -(n+1) * (a/r)^(n +2) * legVec(m+1) * sinVec(m+1);
                hInd = hInd + 1;
            end
        end
    end
    if hInd ~= maxDeg*(maxDeg+2)+1
        error("term matrix is not right.")
    end
    %populate q,s, and Dst dependent terms
    extTerms = [-cosLat , sinLat*cosVec(2), sinLat*sinVec(2);...
        0,sinVec(2),-cosVec(2);...
        sinLat,cosLat*cosVec(2),cosLat*sinVec(2)];
    block(1:3,end-6:end-4) = extTerms;
    %g,q,s Dst-dependent terms
    block(1:3,end-3) = block(1:3,1) * Dst; %g0Dst
    block(1:3,end-2:end) = extTerms * Dst;
end

function F = getIntensity(XYZ)
    %convert vector of B vecs into vector of intensities F
    numMeas = height(XYZ)/3;
    F = zeros(numMeas,1);
    k = 1;
    for n = 1:3:3*numMeas
        F(k) = norm(XYZ(n:n+2));
        k = k+1;
    end
end

function attMatCRF = getMatCRF(q_NEC_CRF)
    %matrix to rotate vector of NEC B vecs into CRF
    A = quat2dcm(q_NEC_CRF);
    A_cell = num2cell(A,[1,2]);
    attMatCRF = blkdiag(A_cell{:});
end

function attMatVFM = getMatVFM(eul,numMeas)
    %matrix to rotate vector of CRF B vecs into VFM
    A_CRF_VFM = angle2dcm(eul(1),eul(2),eul(3),'ZYX');
    C = repmat({A_CRF_VFM},1,numMeas);
    attMatVFM = blkdiag(C{:});
end

function T = getTermMat(Longitude,Latitude,Radius,maxDeg,numModelTerms,Dst)
    %matrix to transform model coefficients into NEC B vecs
    numMeas = length(Longitude);
    w = waitbar(0,"Populating spherical harmonic term matrix.");
    T = zeros(numMeas*3, numModelTerms-3);
    for k = 1:numMeas
        T(3*(k-1)+1:3*(k-1)+3,:) = getTBlock(Longitude(k),Latitude(k),Radius(k),Dst(k),maxDeg,numModelTerms-3);
        waitbar(k/numMeas,w,"Populating spherical harmonic term matrix.");
    end
    close(w);
end

function dataVec = getDataVec(termMat,attMatCRF,modelVec,numMeas)
    %use model vector STRUCT to compute data vector
    eul = modelVec.eul;
    attMatVFM = getMatVFM(eul,numMeas);
    gVec = modelVec.gVec;
    dataXYZ = attMatVFM * attMatCRF * termMat * gVec;
    dataF = getIntensity(dataXYZ);
    dataVec = [dataF;dataXYZ];
end

function dataVec = getDataVec_Vecinput(m,termMat,attMatCRF,numMeas)
    %use model vector AS VECTOR to compute data vector
    eul = m(end-2:end,1);
    gVec = m(1:end-3);
    attMatVFM = getMatVFM(eul,numMeas);
    dataXYZ = attMatVFM * attMatCRF * termMat * gVec;
    dataF = getIntensity(dataXYZ);
    dataVec = [dataF;dataXYZ];
end

function W = getWeightMatrix(Latitude,sigmaScalar,sigmaVector)
    %construct weight matrix for data vector
    cosLat = cosd(Latitude);
    cosLatTrip = [cosLat cosLat cosLat];
    cosLatTrip = cosLatTrip';
    W_s = diag(cosLat/sigmaScalar);
    W_v = diag(cosLatTrip(:)/sigmaVector);
    W = blkdiag(W_s,W_v);
end

function J = parallelJacobian(f, x0)
    f0 = f(x0); 
    M = numel(f0);
    N = numel(x0);
    h = sqrt(eps) * max(abs(x0), 1); 
    J = zeros(M, N);
    parfor(i = 1:N,6)
        x_perturbed = x0; 
        x_perturbed(i) = x_perturbed(i) + h(i);
        f_perturbed = f(x_perturbed);
        J(:, i) = (f_perturbed - f0) / h(i);
    end
end

function [G, err] = approxJacobian(g_i,eul_i,numMeas,termMat,attMatCRF,mode)
    dfunc = @(m) getDataVec_Vecinput(m,termMat,attMatCRF,numMeas);
    if mode == "DERIVEST"
        [G,err] = jacobianest(dfunc,[g_i;eul_i]);
    elseif mode == "FFD"
        G = parallelJacobian(dfunc,[g_i;eul_i]);
        err = [];
    else
        error("Invalid Jacobian approximation mode.")
    end
end