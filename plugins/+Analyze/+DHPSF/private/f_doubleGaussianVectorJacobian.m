% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function jacobian = f_doubleGaussianVectorJacobian(par,ii,jj)
  jacobian1 = exp( -((ii-par(3)).^2+(jj-par(4)).^2) / (2*par(7).^2));
  jacobian2 = exp( -((ii-par(5)).^2+(jj-par(6)).^2) / (2*par(8).^2));
  jacobian3 = par(1) .* exp( -((ii-par(3)).^2+(jj-par(4)).^2) / (2*par(7).^2)) ...
          .* (ii-par(3)) / (par(7).^2);
  jacobian4 = par(1).*exp( -((ii-par(3)).^2+(jj-par(4)).^2) / (2*par(7).^2)) ...
          .* (jj-par(4)) / (par(7).^2);
  jacobian5 = par(2).*exp( -((ii-par(5)).^2+(jj-par(6)).^2) / (2*par(8).^2)) ...
          .* (ii-par(5)) / (par(8).^2);
  jacobian6 = par(2).*exp( -((ii-par(5)).^2+(jj-par(6)).^2) / (2*par(8).^2)) ...
          .* (jj-par(6)) / (par(8).^2);
  jacobian7 = par(1).*exp( -((ii-par(3)).^2+(jj-par(4)).^2) / (2*par(7).^2)) ...
          .* ((ii-par(3)).^2+(jj-par(4)).^2) / (par(7).^3);
  jacobian8 = par(2).*exp( -((ii-par(5)).^2+(jj-par(6)).^2) / (2*par(8).^2)) ...
          .* ((ii-par(5)).^2+(jj-par(6)).^2) / (par(8).^3);
 
  jacobian = [reshape(jacobian1,[],1) reshape(jacobian2,[],1) ...
      reshape(jacobian3,[],1) reshape(jacobian4,[],1) ...
      reshape(jacobian5,[],1) reshape(jacobian6,[],1) ...
      reshape(jacobian7,[],1) reshape(jacobian8,[],1)];
