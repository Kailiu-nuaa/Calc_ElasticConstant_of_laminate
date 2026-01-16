function out = Calc_Laminate_Elastic(EMOD1, EMOD2, ENU12, G12, thetaList, tList)
% Calc_Laminate_Elastic
% Compute in-plane equivalent stiffness of a laminate made of multiple plies.
%
% Inputs:
%   EMOD1, EMOD2 : E1, E2 (Pa)
%   ENU12        : nu12
%   G12          : G12 (Pa)
%   thetaList    : [nPly x 1] ply angles (radians by default)
%   tList        : [nPly x 1] ply thickness (m). Optional: default equal thickness=1
%
% Output struct out:
%   out.C0MAT  : local ply stiffness Q (3x3)
%   out.A      : laminate extensional stiffness (3x3)  (N/m)
%   out.Qeq    : equivalent in-plane stiffness A/h (3x3) (Pa)
%   out.Seq    : equivalent in-plane compliance inv(Qeq) (3x3) (1/Pa)
%   out.Ex, out.Ey, out.Gxy, out.nuxy, out.nuyx
%   out.TRANS_all : 3x3xn transformation matrices
%   out.CMAT_all  : 3x3xn transformed ply stiffness Qbar

    thetaList = thetaList(:);
    nPly = numel(thetaList);

    if nargin < 6 || isempty(tList)
        tList = ones(nPly,1); % default equal thickness = 1 (relative)
    else
        tList = tList(:);
        if numel(tList) ~= nPly
            error('tList must have same length as thetaList.');
        end
    end

    % --- nu21 consistency ---
    ENU21 = ENU12 * EMOD2 / EMOD1;

    % --- local ply stiffness Q (same as your C0MAT) ---
    denom = 1.0 - ENU12 * ENU21;

    C0MAT = zeros(3,3);
    C0MAT(1,1) = EMOD1 / denom;
    C0MAT(2,2) = EMOD2 / denom;
    C0MAT(3,3) = G12;
    C0MAT(1,2) = ENU12 * EMOD2 / denom;
    C0MAT(2,1) = C0MAT(1,2);

    % --- loop plies: build A = sum(Qbar_k * t_k) ---
    A = zeros(3,3);
    TRANS_all = zeros(3,3,nPly);
    CMAT_all  = zeros(3,3,nPly);

    for k = 1:nPly
        TRANS = Transform_Tensor(thetaList(k));
        Qbar  = TRANS * C0MAT * TRANS.';   % same pattern as your CMAT

        TRANS_all(:,:,k) = TRANS;
        CMAT_all(:,:,k)  = Qbar;

        A = A + Qbar * tList(k);
    end

    h = sum(tList);
    Qeq = A / h;        % equivalent in-plane stiffness (Pa)
    Seq = inv(Qeq);     % equivalent compliance (1/Pa)

    % --- engineering constants from Seq ---
    Ex  = 1 / Seq(1,1);
    Ey  = 1 / Seq(2,2);
    Gxy = 1 / Seq(3,3);
    nuxy = -Seq(1,2) / Seq(1,1);
    nuyx = -Seq(1,2) / Seq(2,2);

    % --- pack outputs ---
    out = struct();
    out.C0MAT = C0MAT;
    out.A = A;
    out.h = h;
    out.Qeq = Qeq;
    out.Seq = Seq;
    out.Ex = Ex; out.Ey = Ey; out.Gxy = Gxy;
    out.nuxy = nuxy; out.nuyx = nuyx;
    out.TRANS_all = TRANS_all;
    out.CMAT_all = CMAT_all;
end


function TRANS = Transform_Tensor(theta)
% Same as your Transform_Tensor

    sin_tf = sin(theta);
    cos_tf = cos(theta);

    sin2_tf = sin_tf^2;
    cos2_tf = cos_tf^2;

    TRANS = zeros(3,3);

    TRANS(1,1) = cos2_tf;
    TRANS(2,1) = sin2_tf;
    TRANS(3,1) = sin_tf * cos_tf;

    TRANS(1,2) = sin2_tf;
    TRANS(2,2) = cos2_tf;
    TRANS(3,2) = -sin_tf * cos_tf;

    TRANS(1,3) = -2.0 * sin_tf * cos_tf;
    TRANS(2,3) =  2.0 * sin_tf * cos_tf;
    TRANS(3,3) =  cos2_tf - sin2_tf;
end
