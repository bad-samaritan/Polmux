function varargout=inverse_pmd(brf,options)

%INVERSE_PMD inverse PMD matrix
%   INVERSE_PMD(BRF) applies the inverse PMD matrix of the link to the
%   electric field contained within GSTATE. The link parameter are within
%   the cell BRF. The k-th field of BRF containes the parameters of fiber k
%   of the current link that should be inverted. BRF{k} is a struct with 
%   the following fields:
%
%   BRF{k}.db0 = birefringence [rad] at GSTATE.FN=0 (see RESET_ALL).
%   BRF{k}.theta = azimuth [rad] of all the PMFs composing the fiber.
%   BRF{k}.epsilon = ellepticity [rad] of all the PMFs composing the fiber.
%   BRF{k}.dgdrms = r.m.s. DGD [ns] per trunk.
%   BRF{k}.lcorr = length [m] of each PMF trunk.
%   BRF{k}.betat = beta(omega), i.e. scalar phase shift [rad] including GVD,
%       slope,etc, where omega/2/pi is the vector of FFT frequencies. betat
%       is common to both polarizations.
%   BRF{k}.db1 = differential phase shift [rad] induced by PMD.
%
%   BRF is returned by FIBER.
%
%   UINV=INVERSE_PMD(BRF) also returns in UINV a 3-D matrix such that 
%   UINV(:,:,n) contains  the [2,2] inverse PMD matrix at frequency 
%   GSTATE.FN(n). Such matrix includes also the scalar linear distortion
%   like GVD, with exception (see later).
%
%   [UINV,U]=INVERSE_PMD(BRF) also returns the PMD matrix U.
%
%   INVERSE_PMD(BRF,OPTIONS) allows the following optional parameters:
%
%   OPTIONS.gvd = 'no': do not apply GVD in the U, UINV evaluation, i.e.
%       force BRF{k}.betat = 0.
%   OPTIONS.apply = 'no': do not apply the inverse matrix to GSTATE.FIELDX
%       and GSTATE.FIELDY.
%   OPTIONS.mat = [2,2] unitary matrix that rotates the reference system
%       before applying U. E.g. such matrix can be the one returned by 
%       SET_SOP. 
%
%   Note 1: UINV is a unitary matrix, hence neglects the fiber attenuation.
%           Hence, Uinv(:,:,n) = U(:,:,n)'.
%   Note 2: given the following link, the n-th call to fiber must return 
%           BRF{n}: 
%
%           BRF{1}      BRF{2}      BRF{3}        ...  
%           / \         / \         / \
%           \ /   |\    \ /   |\    \ /   |\      ...   
%   --------------| >---------| >---------| >---------
%                 |/          |/          |/  
%
%   See also FIBER, SET_SOP
%
%   Author: Paolo Serena, 2009
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
%			 
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

global GSTATE

%%%% check
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
if nfc > 1, error('inverse_pmd can be used only with a unique field.'); end

%%%% init
isopt = exist('options','var');
isnotgvd = isopt && isfield(options,'gvd') && strcmp(options.gvd,'no');
Nfft = length(GSTATE.FN);
U = zeros(2,2,Nfft); U(1,1,:) = 1; U(2,2,:) = 1; % PMD matrix for Nfft freq
allgvd = zeros(Nfft,1);
sig0 = eye(2);
sig2 = [0 1;1 0];
sig3i = [0 1;-1 0]; % =i*sig3 = i*[0 -i;i 0]

if isopt && isfield(options,'mat')
    U = update_U(ones(1,Nfft),ones(1,Nfft),options.mat,U);  % change reference system
end

%%%% Calculate matrix U
nfiber = length(brf);
for n=1:nfiber % cycle over the fibers composing the link
    ntrunk = length(brf{n}.theta);
    matR = getmatR(brf{n}.theta(1),brf{n}.epsilon(1),sig0,sig2,sig3i);
    deltabeta=0.5*(brf{n}.db1+brf{n}.db0(1));  % differential beta factor (trunk 1)
    l1 = fastexp(-deltabeta); % eigenvalues
    l2 = 1./l1;
    U = update_U(l1,l2,matR',U);
    for k=2:ntrunk % except the first. The last is finished after.
        matR1 = getmatR(brf{n}.theta(k-1),brf{n}.epsilon(k-1),sig0,sig2,sig3i);
        matR2 = getmatR(brf{n}.theta(k),brf{n}.epsilon(k),sig0,sig2,sig3i);
        matR = matR2'*matR1;

        % Note: Calling A=[GSTATE.FIELDX(k);GSTATE.FIELDY(k)] the electric
        %   field for the kth frequency, we have that matR*D*matR'*A
        %   is the linear PMD step, where D is the diagonal matrix
        %   where the DGD operates.
        % Note 2:
        %
        %   If the overall matrix is R3*D3*R3'*R2*D2*R2'*R1*D1*R1',
        %
        %   l1 and l2 are the eigenvalues of diagonal matrix Dk, 
        %   matR=Rk'*Rk-1, k=2,3.
        
        deltabeta=0.5*(brf{n}.db1+brf{n}.db0(k));  % differential beta factor
        l1 = fastexp(-deltabeta); % eigenvalues
        l2 = 1./l1;
        U = update_U(l1,l2,matR,U);
    end
    matR = getmatR(brf{n}.theta(end),brf{n}.epsilon(end),sig0,sig2,sig3i);
    U = update_U(ones(1,Nfft),ones(1,Nfft),matR,U);  % last trunk
    allgvd = allgvd + brf{n}.betat*brf{n}.lcorr*ntrunk;
end
Hgvd = fastexp(-allgvd);    % overall GVD in one step (scalar component)
if ~isnotgvd
    U(1,1,:) = Hgvd.*squeeze(U(1,1,:)); U(1,2,:) = Hgvd.*squeeze(U(1,2,:)); 
    U(2,1,:) = Hgvd.*squeeze(U(2,1,:)); U(2,2,:) = Hgvd.*squeeze(U(2,2,:)); 
end
Uinv(1,1,:) = conj(U(1,1,:)); Uinv(1,2,:) = conj(U(2,1,:));
Uinv(2,1,:) = conj(U(1,2,:)); Uinv(2,2,:) = conj(U(2,2,:));
if nargout >= 1, varargout{1} = Uinv; end
if nargout >= 2, varargout{2} = U; end

if ~isopt || ~isfield(options,'apply') || strcmp(options.apply,'n')
    uux = fft(GSTATE.FIELDX);
    uuy = fft(GSTATE.FIELDY);
    GSTATE.FIELDX = squeeze(Uinv(1,1,:)).* uux + squeeze(Uinv(1,2,:)).*uuy;
    GSTATE.FIELDY = squeeze(Uinv(2,1,:)).* uux + squeeze(Uinv(2,2,:)).*uuy;
    GSTATE.FIELDX = ifft(GSTATE.FIELDX);
    GSTATE.FIELDY = ifft(GSTATE.FIELDY);
    GSTATE.DISP = zeros(2, GSTATE.NCH);
end   
%--------------------------------------------------------------------------
function U=update_U(l1,l2,matR,Uold)

%UPDATE_U update 3-D matrix U
%   UPDATE_U(L1,L2,MATR,UOLD) updates the 3d matrix UOLD yielding the new
%   one on output U. MATR is the birefringence matrix 2x2 of the new trunk,
%   while L1 and L2 are vectors 1xNfft equal to the eigenvalues. Nfft is
%   the number of frequencies, i.e. the FFT points.
%
%   On output U is o fsize [2,2,Nfft].

matT(1,1,:) = l1*matR(1,1); matT(1,2,:) = l1*matR(1,2);
matT(2,1,:) = l2*matR(2,1); matT(2,2,:) = l2*matR(2,2); % first trunk

U(1,1,:) = matT(1,1,:).*Uold(1,1,:) + matT(1,2,:).*Uold(2,1,:);
U(1,2,:) = matT(1,1,:).*Uold(1,2,:) + matT(1,2,:).*Uold(2,2,:);
U(2,1,:) = -conj(U(1,2,:));%matT(2,1,:).*Uold(1,1,:) + matT(2,2,:).*Uold(2,1,:);
U(2,2,:) = conj(U(1,1,:));%matT(2,1,:).*Uold(1,2,:) + matT(2,2,:).*Uold(2,2,:);
    
%--------------------------------------------------------------------------
function matR=getmatR(theta,epsilon,sig0,sig2,sig3i)

matRth = cos(theta)*sig0 - sin(theta)*sig3i;    % orthogonal matrix
matRepsilon = complex(cos(epsilon)*sig0,sin(epsilon)*sig2); % unitary
matR = matRth*matRepsilon;  % matrix of change of basis over the PSPs.
