Browse : [JavaScript](https://github.com/michel-leonard/ciede2000-javascript) · [Java](https://github.com/michel-leonard/ciede2000-java) · [Julia](https://github.com/michel-leonard/ciede2000-julia) · [Kotlin](https://github.com/michel-leonard/ciede2000-kotlin) · [Lua](https://github.com/michel-leonard/ciede2000-lua) · **MATLAB** · [Microsoft Excel](https://github.com/michel-leonard/ciede2000-excel) · [PHP](https://github.com/michel-leonard/ciede2000-php) · [Perl](https://github.com/michel-leonard/ciede2000-perl) · [Python](https://github.com/michel-leonard/ciede2000-python) · [R](https://github.com/michel-leonard/ciede2000-r)

# CIEDE2000 color difference formula in MATLAB

This page presents the CIEDE2000 color difference, implemented in the MATLAB programming language.

![Logo of CIEDE2000 in MATrix LABoratory](https://raw.githubusercontent.com/michel-leonard/ciede2000-color-matching/refs/heads/main/docs/assets/images/logo.jpg)

## About

Here you’ll find the first rigorously correct implementation of CIEDE2000 that doesn’t use any conversion between degrees and radians. Set parameter `canonical` to obtain results in line with your existing pipeline.

`canonical`|The algorithm operates...|
|:--:|-|
`false`|in accordance with the CIEDE2000 values currently used by many industry players|
`true`|in accordance with the CIEDE2000 values provided by [this](https://hajim.rochester.edu/ece/sites/gsharma/ciede2000/) academic MATLAB function|

## Our CIEDE2000 offer

This production-ready file, released in 2026, contain the CIEDE2000 algorithm.

Source File|Type|Bits|Purpose|Advantage|
|:--:|:--:|:--:|:--:|:--:|
[ciede2000.m](./ciede2000.m)|`double`|64|General|Vectorization, Interoperability|

`ciede2000` has a low memory cost, allocating no more than 9 temporary vectors at any one time.

### Software Versions

- MATLAB 25.2
- GNU Octave 8.4

### Example Usage

We calculate the CIEDE2000 distance between two colors, first without and then with parametric factors.

```matlab
% You must name this file "example.m" and place it in the same directory as "ciede2000.m"

function example()
	% Example of two L*a*b* colors
	l1 = 43.8; a1 = 57.1; b1 = -5.2;
	l2 = 44.1; a2 = 91.2; b2 = 13.8;

	delta_e = ciede2000(l1, a1, b1, l2, a2, b2);
	fprintf("delta_e = %.14f\n", delta_e);
	% ΔE2000 = 10.71979501772133

	% Example of parametric factors used in the dental industry
	kl = 2.0;  kc = 1.0;  kh = 1.0;

	% Perform a CIEDE2000 calculation compliant with that of Gaurav Sharma
	canonical = true;

	delta_e = ciede2000(l1, a1, b1, l2, a2, b2, kl, kc, kh, canonical);
	fprintf("delta_e = %.14f\n", delta_e);
	% ΔE2000 = 10.71709331390274
end
```

**Note**: this example uses scalars, but if you prefer, it would have worked just as well with vectors.

You can then run `octave-cli example.m` to display the calculated CIEDE2000 color differences.

### Test Results

LEONARD’s tests are based on well-chosen L\*a\*b\* colors, with various parametric factors `kL`, `kC` and `kH`.

```
CIEDE2000 Verification Summary :
          Compliance : [ ] CANONICAL [X] SIMPLIFIED
  First Checked Line : 30,32,32,30,128,-127.40000000000001,1,1,1,45.320587416812415
           Precision : 12 decimal digits
           Successes : 100000000
               Error : 0
            Duration : 301.53 seconds
     Average Delta E : 67.13
   Average Deviation : 6.4e-15
   Maximum Deviation : 3.1e-13
```

```
CIEDE2000 Verification Summary :
          Compliance : [X] CANONICAL [ ] SIMPLIFIED
  First Checked Line : 30,32,32,30,128,-127.40000000000001,1,1,1,45.320397256377987
           Precision : 12 decimal digits
           Successes : 100000000
               Error : 0
            Duration : 300.26 seconds
     Average Delta E : 67.13
   Average Deviation : 6.7e-15
   Maximum Deviation : 3.1e-13
```

## Public Domain Licence

You are free to use these files, even for commercial purposes.
