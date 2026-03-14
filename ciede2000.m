% This function written in MATLAB is not affiliated with the CIE (International Commission on Illumination),
% and is released into the public domain. It is provided "as is" without any warranty, express or implied.

function delta_e = ciede2000_classic(l1, a1, b1, l2, a2, b2, kl, kc, kh, canonical)
	% Sets the default value for the last 4 parameters, which are optional.
	if nargin < 9 kl = 1.0; kc = 1.0; kh = 1.0; end;  % scalars.
	canonical = 9 < nargin && canonical;              % boolean.
	% This scalar expansion wrapper works with numbers, not vectors.
	delta_e = ciede2000([l1], [a1], [b1], [l2], [a2], [b2], [kl], [kc], [kh], canonical);
	delta_e = delta_e(1);
end

% The classic vectorized CIE ΔE2000 implementation, which operates on two L*a*b* colors, and returns their difference.
% "l" ranges from 0 to 100, while "a" and "b" are unbounded and commonly clamped to the range of -128 to 127.
function delta_e = ciede2000(l1, a1, b1, l2, a2, b2, kl, kc, kh, canonical)
	% Sets the default value for the last 4 parameters, which are optional.
	if nargin < 9 kl = 1.0; kc = 1.0; kh = 1.0; end;  % scalars or vectors.
	canonical = 9 < nargin && canonical;              % boolean.
	% Working in MATLAB with the CIEDE2000 color-difference formula.
	% kl, kc, kh are parametric factors to be adjusted according to
	% different viewing parameters such as textures, backgrounds...
	n = ((sqrt(a1 .* a1 + b1 .* b1) + sqrt(a2 .* a2 + b2 .* b2)) * 0.5) .^ 7.0;
	% A factor involving chroma raised to the power of 7 designed to make
	% the influence of chroma on the total color difference more accurate.
	n = 1.0 + 0.5 * (1.0 - sqrt(n ./ (n + 6103515625.0)));
	% Application of the chroma correction factor.
	c1 = sqrt(a1 .* a1 .* n .* n + b1 .* b1);
	c2 = sqrt(a2 .* a2 .* n .* n + b2 .* b2);
	% atan2 is preferred over atan because it accurately computes the angle of
	% a point (x, y) in all quadrants, handling the signs of both coordinates.
	h1 = atan2(b1, a1 .* n);
	h2 = atan2(b2, a2 .* n);
	% Vectorized conditionals
	mask = h1 < 0.0;
	h1(mask) = h1(mask) + 2.0 * pi;
	mask = h2 < 0.0;
	h2(mask) = h2(mask) + 2.0 * pi;
	% When the hue angles lie in different quadrants, the straightforward
	% average can produce a mean that incorrectly suggests a hue angle in
	% the wrong quadrant, the next lines handle this issue.
	h_mean = (h1 + h2) * 0.5;
	h_delta = (h2 - h1) * 0.5;
	% Vectorized conditionals
	mask = pi + 1E-14 < abs(h2 - h1);
	h1 = false; h2 = false; % Not used anymore.
	h_delta(mask) = h_delta(mask) + pi;
	% canonical = 0 ==> Lindbloom’s implementation, Netflix’s VMAF, ...
	% canonical = 1 ==> Sharma’s implementation, OpenJDK, ...
	n = canonical & (pi + 1E-14 < h_mean(mask));
	h_mean(mask) = h_mean(mask) + pi * (~n - n);
	p = 36.0 * h_mean - 55.0 * pi;
	n = ((c1 + c2) * 0.5) .^ 7.0;
	% The hue rotation correction term is designed to account for the
	% non-linear behavior of hue differences in the blue region.
	r_t = -2.0 * sqrt(n ./ (n + 6103515625.0));
	r_t = r_t .* sin(pi / 3.0 * exp((p .* p) / (-25.0 * pi * pi)));
	mask = false; p = false; % Not used anymore.
	n = ((l1 + l2) * 0.5 - 50.0) .^ 2.0;
	% Lightness.
	l = (l2 - l1) ./ (kl .* (1.0 + 0.015 * n ./ sqrt(20.0 + n)));
	% These coefficients adjust the impact of different harmonic
	% components on the hue difference calculation.
	t = 1.0 - 0.17 .* sin(h_mean + pi / 3.0);
	t = t + 0.24 * sin(2.0 * h_mean + pi * 0.5);
	t = t + 0.32 * sin(3.0 * h_mean + 8.0 * pi / 15.0);
	t = t - 0.20 * sin(4.0 * h_mean + 3.0 * pi / 20.0);
	n = c1 + c2;
	% Hue.
	h = 2.0 * sqrt(c1 .* c2) .* sin(h_delta);
	h = h ./ (kh .* (1.0 + 0.0075 * n .* t));
	h_mean = false; h_delta = false; t = false; % Not used anymore.
	% Chroma.
	c = (c2 - c1) ./ (kc .* (1.0 + 0.0225 * n));
	% The result reflects the actual geometric distance in the color space, given a tolerance of 1.3e-12.
	delta_e = sqrt(l .* l + h .* h + c .* c + c .* h .* r_t);
end

% If you remove the constant 1E-14, the code will continue to work, but CIEDE2000
% interoperability between all programming languages will no longer be guaranteed.
