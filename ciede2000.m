% This function written in MATLAB is not affiliated with the CIE (International Commission on Illumination),
% and is released into the public domain. It is provided "as is" without any warranty, express or implied.

% The classic vectorized CIE ΔE2000 implementation, which operates on two L*a*b* colors, and returns their difference.
% "l" ranges from 0 to 100, while "a" and "b" are unbounded and commonly clamped to the range of -128 to 127.
function delta_e = ciede2000(l1, a1, b1, l2, a2, b2, kl, kc, kh, canonical)
	% Sets the default value for the last 4 parameters, which are optional.
	if nargin < 9 kl = 1.0; kc = 1.0; kh = 1.0; end;  % scalars or vectors.
	canonical = 9 < nargin && canonical;              % boolean.
	% Working in MATLAB with the CIEDE2000 color-difference formula.
	% kl, kc, kh are parametric factors to be adjusted according to
	% different viewing parameters such as textures, backgrounds...
	c1 = b1 .* b1;
	c2 = b2 .* b2;
	n = sqrt(a1 .* a1 + c1);
	n = n + sqrt(a2 .* a2 + c2);
	n = n * 0.5;
	n = n .^ 7.0;
	n = n ./ (n + 6103515625.0);
	n = sqrt(n);
	n = n * 0.5;
	n = 1.5 - n; % The chroma correction factor.

	% atan2 is preferred over atan because it accurately computes the angle of
	% a point (x, y) in all quadrants, handling the signs of both coordinates.
	c = a1 .* n;
	c1 = c1 + c .* c;
	c1 = sqrt(c1);
	hm = atan2(b1, c);

	c = a2 .* n;
	c2 = c2 + c .* c;
	c2 = sqrt(c2);
	h = atan2(b2, c);

	n = hm < 0.0;
	hm(n) = hm(n) + 2.0 * pi;

	n = h < 0.0;
	h(n) = h(n) + 2.0 * pi;

	% When the hue angles lie in different quadrants, the straightforward
	% average can produce a mean that incorrectly suggests a hue angle in
	% the wrong quadrant, the next 10 lines handle this issue.
	h = h - hm;
	n = pi + 1E-14 < abs(h);
	h = h * 0.5; % h_delta
	hm = hm + h; % h_mean

	% The part where most programmers get it wrong.
	h(n) = h(n) + pi;
	% canonical = 0 ==> Lindbloom’s implementation, Netflix’s VMAF, ...
	% canonical = 1 ==> Sharma’s implementation, OpenJDK, ...
	c = canonical & (pi + 1E-14 < hm(n));
	hm(n) = hm(n) + pi * (~c - c);

	% The hue rotation correction term is designed to account for the
	% non-linear behavior of hue differences in the blue region.
	c = c1 + c2;
	n = c * 0.5;
	n = n .^ 7.0;
	n = n ./ (n + 6103515625.0);
	n = sqrt(n);
	r_t = -2.0 * n;
	n = 36.0 * hm;
	n = n - 55.0 * pi;
	n = n / (5.0 * pi);
	n = n .* n;
	n = exp(-n);
	n = n * (pi / 3.0);
	n = sin(n);
	r_t = r_t .* n;

	% These coefficients adjust the impact of different harmonic
	% components on the hue difference calculation.
	n = 150.0;
	n = n - 25.5 * sin(hm + pi / 3.0);
	n = n + 36.0 * sin(2.0 * hm + pi * 0.5);
	n = n + 48.0 * sin(3.0 * hm + 8.0 * pi / 15.0);
	n = n - 30.0 * sin(4.0 * hm + 3.0 * pi / 20.0);
	n = n / 20000.0;
	n = n .* c;
	n = n + 1.0;
	n = n .* kh;
	% Hue.
	h = 2.0 * sin(h);
	h = h .* sqrt(c1 .* c2);
	h = h ./ n;
	hm = false; % Not used anymore.

	% Lightness.
	l = l2 - l1;
	n = l1 + l2;
	n = n * 0.5;
	n = n - 50.0;
	n = n .* n;
	n = n ./ sqrt(20.0 + n);
	n = n * 3.0;
	n = n / 200.0;
	n = n + 1.0;
	n = n .* kl;
	l = l ./ n;

	% Chroma.
	c = c * 9.0;
	c = c / 400.0;
	c = c + 1.0;
	c = c .* kc;
	c = (c2 - c1) ./ c;
	n = false; c1 = false; c2 = false; % Not used anymore.

	delta_e = l .* l;
	delta_e = delta_e + h .* h;
	delta_e = delta_e + c .* c;
	h = h .* r_t;
	delta_e = delta_e + c .* h;
	% The result reflects the actual geometric distance in the color space, given a tolerance of 1.3e-12.
	% The function allocates no more than 9 temporary vectors at any one time.
	delta_e = sqrt(delta_e);
end

% If you remove the constant 1E-14, the code will continue to work, but CIEDE2000
% interoperability between all programming languages will no longer be guaranteed.

% Source code tested by Michel LEONARD
% Website: ciede2000.pages-perso.free.fr
