clear
reset
set key off
set border 3

# Add a vertical dotted line at x=0 to show centre (mean) of distribution.
set yzeroaxis

# Each bar is half the (visual) width of its x-range.

bin_width = 100;
bin_number(x) = floor(x/bin_width)

set boxwidth 0.9 * bin_width
set style fill solid 0.5

rounded(x) = bin_width * (bin_number(x) + 0.5)

plot '../bin/path-count.report' using (rounded($4)):(1) smooth frequency with boxes

