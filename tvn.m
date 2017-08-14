function p = tvn( lw, up, cr, epsi ) 
%TVN
%  A function for computing trivariate normal probabilities.
%    p = tvn( lw, up, cr, epsi ) 
%  Parameters
%     lw  array of lower integration limits.
%     up  array of upper integration limits.
%     cr   array of correlation coefficients, in order c21, c31, c32.
%     epsi  optional (default = 1e-7) absolute accuracy; highest accuracy 
%            for most computations is approximately 1e-14.
%  Example 
%   p = tvn( [-3,-4,-inf], [inf,4,3], [.3,-.4,.5] ) 
%

%  This file contains functions tvn, tvnl (trivariate normal), bvn,
%  bvnl (bivariate normal), phid (univariate normal), plus support functions.
%
%   Copyright (C) 2011, Alan Genz,  All rights reserved.               
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided the following conditions are met:
%     1. Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%     2. Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%     3. The contributor name(s) may not be used to endorse or promote 
%        products derived from this software without specific prior written 
%        permission.
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
%   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
%   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
  if nargin < 4, epsi = 1e-7; end, ep = max( [ 1e-14 epsi ] ); 
  [m n] = size(lw); if m==1, lu = [lw; up]; else, lu = [lw'; up']; end  
  tvn = tvnl( up, cr, ep ) -  tvnl( lw, cr, ep ); 
  tvn = tvn - tvnl( [lu(1,1) lu(2,2) lu(2,3)], cr, ep ); 
  tvn = tvn - tvnl( [lu(2,1) lu(1,2) lu(2,3)], cr, ep ); 
  tvn = tvn - tvnl( [lu(2,1) lu(2,2) lu(1,3)], cr, ep ); 
  tvn = tvn + tvnl( [lu(1,1) lu(1,2) lu(2,3)], cr, ep ); 
  tvn = tvn + tvnl( [lu(1,1) lu(2,2) lu(1,3)], cr, ep ); 
  tvn = tvn + tvnl( [lu(2,1) lu(1,2) lu(1,3)], cr, ep ); 
  p = min( max( 0, tvn ), 1 );
%
function p = tvnl( h, r, epsi )
%TVNL    
%     A function for computing trivariate normal probabilities.
%     It calculates the probability that x(i) < h(i), for i = 1, 2, 3.     
%    h     real array of three upper limits for probability distribution 
%    r     real array of three correlation coefficients, r should 
%          contain the lower left portion of the correlation matrix R. 
%          r should contain the values r21, r31, r32 in that order.
%    epsi  optional (default = 1e-7) absolute accuracy; maximum accuracy 
%            for most computations is approximately 1e-14.
%   Example: p = tvnl( [1 2 4], [.3 .4 -.6] )
% 

%     This function uses algorithms developed from the ideas 
%     described in the papers:
%       R.L. Plackett, Biometrika 41(1954), pp. 351-360.
%       Z. Drezner, Math. Comp. 62(1994), pp. 289-294.
%     with adaptive integration from (0,0,r32) for R. 
%     The software is based on work described in the paper
%      "Numerical Computation of Rectangular Bivariate and Trivariate 
%        Normal and t Probabilities", by
%          Alan Genz
%          Department of Mathematics
%          Washington State University
%          Pullman, WA 99164-3113
%          Email : alangenz@wsu.edu
%     This file contains functions tvnl (trivariate normal), bvnl 
%      (bivariate normal), phid (univariate normal), plus support functions.
%
%   Copyright (C) 2011, Alan Genz,  All rights reserved.               
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided the following conditions are met:
%     1. Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%     2. Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%     3. The contributor name(s) may not be used to endorse or promote 
%        products derived from this software without specific prior written 
%        permission.
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
%   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
%   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
  if nargin < 3, epsi = 1e-7; end, epst = max( [ 1e-14 epsi ] ); 
  if h(1) == -inf | h(2) == -inf | h(3) == -inf, p = 0;
  elseif h(1) == inf
    if h(2) == inf, if h(3) == inf, p = 1; else p = phid(h(3)); end
    elseif h(3) == inf, p = phid(h(2)); 
    else p = bvnl( h(2), h(3), r(3) ); 
    end
  elseif h(2) == inf
    if h(3) == inf, p = phid(h(1)); 
    else p = bvnl( h(1), h(3), r(2) ); 
    end
  elseif h(3) == inf, p = bvnl( h(1), h(2), r(1) );
  else h1 = h(1); h2 = h(2); h3 = h(3); r12 = r(1); r13 = r(2); r23 = r(3);
    %
    %     Sort R's and check for special cases
    %
    if abs(r12) > abs(r13), h2 = h3; h3 = h(2); r12 = r13; r13 = r(1); end
    if abs(r13) > abs(r23), h1 = h2; h2 = h(1); r23 = r13; r13 = r(3); end
    if abs(h1) + abs(h2) + abs(h3) < epst 
      tvn = ( 1 + 2*( asin(r12) + asin(r13) + asin(r23) )/pi )/8;
    elseif abs(r12) + abs(r13) < epst, tvn = phid(h1)*bvnl( h2, h3, r23 );
    elseif abs(r13) + abs(r23) < epst, tvn = phid(h3)*bvnl( h1, h2, r12 );
    elseif abs(r12) + abs(r23) < epst, tvn = phid(h2)*bvnl( h1, h3, r13 );
    elseif 1 - r23 < epst, tvn = bvnl( h1, min( [h2 h3] ), r12 );
    elseif r23 + 1 < epst & h2 > -h3 
      tvn = bvnl( h1, h2, r12 ) - bvnl( h1, -h3, r12 );
    else a12 = asin(r12); a13 = asin(r13);
      %
      %        Use numerical integration to compute probability
      %           and add singular value
      %
      tvn = adonet( @(x)tvnf( x, h1,h2,h3,r23,a12,a13 ), 0, 1, epst )/(2*pi);
      tvn = tvn + bvnl( h2, h3, r23 )*phid(h1); 
    end
    p = max( 0, min( tvn, 1 ) );
  end
%
% end tvnl
%
function f = tvnf( x,  h1, h2, h3, r23, a12, a13 )
%
%     Computes Plackett formula integrands
%
  f = 0; [ r12 rr2 ] = sincs( a12*x ); [ r13 rr3 ] = sincs( a13*x );
  if abs(a12) > 0, f = f + a12*pntgnd( h1, h2, h3, r13, r23, r12, rr2 ); end
  if abs(a13) > 0, f = f + a13*pntgnd( h1, h3, h2, r12, r23, r13, rr3 ); end
%
% end tvnf
%
function [sx, cs] = sincs( x ), ee = ( pi/2 - abs(x) )^2;
%
%     Computes sin(x), cos(x)^2, with series approx. for |x| near pi/2
%
  if ee < 5e-5, cs = ee*( 1 - ee*( 1 - 2*ee/15 )/3 );
    sx = ( 1 - ee*( 1 - ee/12 )/2 )*sign( x );
  else, sx = sin(x); cs = 1 - sx*sx;
  end
%
% end sincs
%
function f = pntgnd( ba, bb, bc, ra, rb, r, rr )
%
%     Computes Plackett formula integrand;
%
  f = 0; dt = rr*( rr - ( ra - rb )^2 - 2*ra*rb*( 1 - r ) );
  if dt > 0  
    bt = ( bc*rr + ba*( r*rb - ra ) + bb*( r*ra -rb ) )/sqrt(dt) ;
    ft = ( ba - r*bb )^2/rr + bb*bb;
    if bt > -10 & ft < 100, f = exp( -ft/2 );
      if bt < 10,  f = f*phid(bt); end
    end
  end
%
% end pntgnd
%
function p = phid(z), p = erfc(-z/sqrt(2))/2; % Normal cdf
%
% end phid
%
function p = bvnl( dh, dk, r )
%BVNL
  p = bvnu( -dh, -dk, r );
%
%   end bvnl
%
function p = bvnu( dh, dk, r )
%
%BVNU
%  A function for computing bivariate normal probabilities.
%  bvnu calculates the probability that x > dh and y > dk. 
%    parameters  
%      dh 1st lower integration limit
%      dk 2nd lower integration limit
%      r   correlation coefficient
%  Example: p = bvnu( -3, -1, .35 )
%  Note: to compute the probability that x < dh and y < dk, 
%        use bvnu( -dh, -dk, r ). 
%

%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
%    This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for Matlab by Alan Genz. Minor bug modifications 7/98, 2/10.
%
  if dh ==  inf | dk ==  inf, p = 0;
  elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
  elseif dk == -inf, p = phid(-dh);
  elseif r == 0, p = phid(-dh)*phid(-dk);  
  else, tp = 2*pi; h = dh; k = dk; hk = h*k; bvn = 0; 
    if abs(r) < 0.3      % Gauss Legendre points and weights, n =  6
      w(1:3) = [0.1713244923791705 0.3607615730481384 0.4679139345726904];
      x(1:3) = [0.9324695142031522 0.6612093864662647 0.2386191860831970];
    elseif abs(r) < 0.75 % Gauss Legendre points and weights, n = 12
      w(1:3) = [.04717533638651177 0.1069393259953183 0.1600783285433464];
      w(4:6) = [0.2031674267230659 0.2334925365383547 0.2491470458134029];
      x(1:3) = [0.9815606342467191 0.9041172563704750 0.7699026741943050];
      x(4:6) = [0.5873179542866171 0.3678314989981802 0.1252334085114692];
    else,                % Gauss Legendre points and weights, n = 20
      w(1:3) = [.01761400713915212 .04060142980038694 .06267204833410906];
      w(4:6) = [.08327674157670475 0.1019301198172404 0.1181945319615184];
      w(7:9) = [0.1316886384491766 0.1420961093183821 0.1491729864726037];
      w(10) =   0.1527533871307259;
      x(1:3) = [0.9931285991850949 0.9639719272779138 0.9122344282513259];
      x(4:6) = [0.8391169718222188 0.7463319064601508 0.6360536807265150];
      x(7:9) = [0.5108670019508271 0.3737060887154196 0.2277858511416451];
      x(10) =   0.07652652113349733;
    end, w = [w  w]; x = [1-x 1+x]; 
    if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r)/2;  
      sn = sin(asr*x); bvn = exp((sn*hk-hs)./(1-sn.^2))*w';
      bvn = bvn*asr/tp + phid(-h)*phid(-k);  
    else, if r < 0, k = -k; hk = -hk; end
      if abs(r) < 1, as = 1-r^2; a = sqrt(as); bs = (h-k)^2;
        asr = -( bs/as + hk )/2; c = (4-hk)/8 ; d = (12-hk)/80; 
        if asr > -100, bvn = a*exp(asr)*(1-c*(bs-as)*(1-d*bs)/3+c*d*as^2); end
        if hk  > -100, b = sqrt(bs); sp = sqrt(tp)*phid(-b/a);
          bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1-d*bs)/3 );
        end, a = a/2; xs = (a*x).^2; asr = -( bs./xs + hk )/2; 
        ix = find( asr > -100 ); xs = xs(ix); sp = ( 1 + c*xs.*(1+5*d*xs) ); 
        rs = sqrt(1-xs); ep = exp( -(hk/2)*xs./(1+rs).^2 )./rs; 
        bvn = ( a*( (exp(asr(ix)).*(sp-ep))*w(ix)' ) - bvn )/tp; 
      end 
      if r > 0, bvn =  bvn + phid( -max( h, k ) ); 
      elseif h >= k, bvn = -bvn;
      else, if h < 0, L = phid(k)-phid(h); else, L = phid(-h)-phid(-k); end
        bvn =  L - bvn;
      end
    end, p = max( 0, min( 1, bvn ) );
  end
%
%   end bvnu
%
function fin = adonet( f, a, b, tol ), nl = 100;
%
%     1-dimensional adaptive integration
%
  ai(1) = a; bi(1) = b; ip = 1; im = 1; err = 1;
  while 4*err > tol & im < nl, im = im + 1;
    bi(im) = bi(ip); ai(im) = ( ai(ip) + bi(ip) )/2; bi(ip) = ai(im);
    [ fi(ip) ei(ip) ] = krnrdt( ai(ip), bi(ip), f );
    [ fi(im) ei(im) ] = krnrdt( ai(im), bi(im), f );
    fin = sum( fi(1:im) ); err = sqrt( sum( ei(1:im).^2 ) ); 
    [ ex, ip ] = max(ei);
  end
%
function [ resk, err ] = krnrdt( a, b, f )
%
%     Kronrod integration rule
%
  wg0    = 0.2729250867779007;
  wg(1:3) = [ .05566856711617449 0.1255803694649048 0.1862902109277352 ];
  wg(4:5) = [ 0.2331937645919914 0.2628045445102478 ];
%
  xgk(1:3) = [ 0.9963696138895427 0.9782286581460570 0.9416771085780681 ];
  xgk(4:6) = [ 0.8870625997680953 0.8160574566562211 0.7301520055740492 ];
  xgk(7:9) = [ 0.6305995201619651 0.5190961292068118 0.3979441409523776 ];
  xgk(10:11) = [ 0.2695431559523450 0.1361130007993617 ];
%
  wgk0 = 0.1365777947111183;
  wgk(1:3) = [ .00976544104596129 .02715655468210443 .04582937856442671 ];
  wgk(4:6) = [ .06309742475037484 .07866457193222764 .09295309859690074 ];
  wgk(7:9) = [ 0.1058720744813894 0.1167395024610472 0.1251587991003195 ];
  wgk(10:11) = [ 0.1312806842298057 0.1351935727998845 ];
  wid = ( b - a )/2; cen = ( b + a )/2; fc = feval( f, cen );
  resg = fc*wg0; resk = fc*wgk0;
  for j = 1 : 5, t = wid*xgk(2*j-1); 
    fc = feval( f, cen - t ) + feval( f, cen + t );
    resk = resk + wgk(2*j-1)*fc; t = wid*xgk(2*j); 
    fc = feval( f, cen - t ) + feval( f, cen + t );
    resk = resk + wgk(2*j)*fc; resg = resg + wg(j)*fc; 
  end
  t = wid*xgk(11); fc = feval( f, cen - t ) + feval( f, cen + t );
  resk = wid*( resk + wgk(11)*fc ); err = abs( resk - wid*resg  );
%
%     end krnrdt
%